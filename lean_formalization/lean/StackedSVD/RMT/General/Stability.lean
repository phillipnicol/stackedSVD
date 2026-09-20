/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.R4C
import StackedSVD.RMT.MP7

/-!
# Item G3a: deterministic stability lemmas for the general-law trace step

`notes/archive/prop_single_table_general.md` section 5, unit G3a. This file holds the deterministic
stability lemmas that unit G3 (`RMT/General/Trace.lean`) needs for its self-consistent
equation step. Every declaration below is a copy of a declaration of `RMT/R1.lean`, under the
same name, moved into the `GenRMT` namespace. Each docstring names the source file and line.

The copy exists because choice 8 of the plan note keeps every file under `RMT/General/` free
of `RMT/R1.lean`, `RMT/R2.lean`, `RMT/SteinStep.lean` and `Vendor/COLT83/`. Gate 8
(`scripts/check_core_imports.py`) checks that boundary. All twelve declarations are
deterministic: they hold for every real symmetric matrix `W` and carry no probability and no
Gaussian argument. Item F34 of `notes/FOLLOWUP_LIST.md` tracks the dedup of this copy against
its `R1` original.
-/

open Filter Topology Set
open scoped Matrix

namespace StackedSVD
namespace GenRMT

variable {d : ℕ} {W : Matrix (Fin d) (Fin d) ℝ} {z : ℂ} {c : ℝ}

/-- The `δ` of step 4, a lower bound on `‖E s_N‖`.

Copy of `R1.rootLb` (`RMT/R1.lean:162`); choice 8 keeps `RMT/General/` free of that file,
item F34. -/
noncomputable def rootLb (c : ℝ) (z : ℂ) : ℝ := min 1 (1 / (2 * (‖z‖ + ‖z + 1 - c‖)))

/-- Copy of `R1.rootLb_pos` (`RMT/R1.lean:446`); choice 8 keeps `RMT/General/` free of that
file, item F34. -/
theorem rootLb_pos (hz : z ≠ 0) : 0 < rootLb c z := by
  have hzn : 0 < ‖z‖ := norm_pos_iff.mpr hz
  have hnn : (0 : ℝ) ≤ ‖z + 1 - (c : ℂ)‖ := norm_nonneg _
  refine lt_min one_pos ?_
  positivity

/-- `Im s ≥ η ‖s‖²`, by Cauchy-Schwarz on the eigenvalue sum.

Copy of `R1.im_stieltjesC_ge` (`RMT/R1.lean:183`); choice 8 keeps `RMT/General/` free of that
file, item F34. -/
theorem im_stieltjesC_ge (hW : W.IsHermitian) (hz : 0 < z.im) :
    z.im * ‖R4C.stieltjesC W z‖ ^ 2 ≤ (R4C.stieltjesC W z).im := by
  rcases Nat.eq_zero_or_pos d with hd0 | hd0
  · subst hd0
    simp [R4C.stieltjesC]
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd0
  set r : Fin d → ℝ := fun a => ‖((hW.eigenvalues a : ℂ) - z)⁻¹‖ with hrdef
  have hterm : ∀ a : Fin d, (((hW.eigenvalues a : ℂ) - z)⁻¹).im = z.im * r a ^ 2 := by
    intro a
    have him : ((hW.eigenvalues a : ℂ) - z).im = -z.im := by simp
    have hra : r a ^ 2 = (Complex.normSq ((hW.eigenvalues a : ℂ) - z))⁻¹ := by
      simp only [hrdef, norm_inv, inv_pow, Complex.sq_norm]
    rw [Complex.inv_im, him, hra]
    ring
  have hScb : ‖(R4C.resolvC W z).trace‖ ≤ ∑ a, r a := by
    rw [R4C.trace_resolvC hW hz.ne']
    exact norm_sum_le _ _
  have hnd : ‖(d : ℂ)‖ = (d : ℝ) := by simp
  have hnorm : ‖R4C.stieltjesC W z‖ ≤ (d : ℝ)⁻¹ * ∑ a, r a := by
    rw [R4C.stieltjesC, norm_mul, norm_inv, hnd]
    exact mul_le_mul_of_nonneg_left hScb (by positivity)
  have him : (R4C.stieltjesC W z).im = (d : ℝ)⁻¹ * (z.im * ∑ a, r a ^ 2) := by
    rw [R4C.stieltjesC, R4C.trace_resolvC hW hz.ne']
    have hdc : ((d : ℕ) : ℂ)⁻¹ = (((d : ℝ)⁻¹ : ℝ) : ℂ) := by push_cast; ring
    rw [hdc]
    simp only [Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im, zero_mul, add_zero]
    congr 1
    rw [Complex.im_sum, Finset.mul_sum]
    exact Finset.sum_congr rfl fun a _ => hterm a
  have hcs : (∑ a, r a) ^ 2 ≤ (d : ℝ) * ∑ a, r a ^ 2 := by
    have h := sq_sum_le_card_mul_sum_sq (s := (Finset.univ : Finset (Fin d))) (f := r)
    simpa using h
  have hrnn : 0 ≤ ∑ a, r a := Finset.sum_nonneg fun a _ => norm_nonneg _
  rw [him]
  calc z.im * ‖R4C.stieltjesC W z‖ ^ 2
      ≤ z.im * ((d : ℝ)⁻¹ * ∑ a, r a) ^ 2 :=
        mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (norm_nonneg _) hnorm 2) hz.le
    _ = z.im * ((d : ℝ)⁻¹) ^ 2 * (∑ a, r a) ^ 2 := by ring
    _ ≤ z.im * ((d : ℝ)⁻¹) ^ 2 * ((d : ℝ) * ∑ a, r a ^ 2) :=
        mul_le_mul_of_nonneg_left hcs (by positivity)
    _ = (d : ℝ)⁻¹ * (z.im * ∑ a, r a ^ 2) := by
        field_simp

/-- `d⁻¹ tr G(z)²` is the derivative in `z` of `d⁻¹ tr G(z)`.

Copy of `R1.hasDerivAt_stieltjesC` (`RMT/R1.lean:227`); choice 8 keeps `RMT/General/` free of
that file, item F34. -/
theorem hasDerivAt_stieltjesC (hW : W.IsHermitian) (hz : 0 < z.im) :
    HasDerivAt (R4C.stieltjesC W) (R4C.stieltjes2C W z) z := by
  have hopen : IsOpen {ζ : ℂ | 0 < ζ.im} := isOpen_lt continuous_const Complex.continuous_im
  have hmem : z ∈ {ζ : ℂ | 0 < ζ.im} := hz
  have hkey : HasDerivAt (fun ζ : ℂ => (d : ℂ)⁻¹ * ∑ a, ((hW.eigenvalues a : ℂ) - ζ)⁻¹)
      ((d : ℂ)⁻¹ * ∑ a, (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2) z := by
    refine HasDerivAt.const_mul _ (HasDerivAt.fun_sum fun a _ => ?_)
    have h1 : HasDerivAt (fun ζ : ℂ => (hW.eigenvalues a : ℂ) - ζ) (-1) z := by
      have h0 : HasDerivAt (fun ζ : ℂ => (hW.eigenvalues a : ℂ) - ζ) (-(1 : ℂ)) z :=
        HasDerivAt.const_sub ((hW.eigenvalues a : ℂ)) (hasDerivAt_id z)
      exact h0
    have h2 := h1.inv (R4C.eigenvalue_sub_ne_zero hW hz.ne' a)
    have heq : -(-1 : ℂ) / ((hW.eigenvalues a : ℂ) - z) ^ 2
        = (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2 := by
      rw [inv_pow, neg_neg, one_div]
    rw [← heq]
    exact h2
  have hval : R4C.stieltjes2C W z = (d : ℂ)⁻¹ * ∑ a, (((hW.eigenvalues a : ℂ) - z)⁻¹) ^ 2 := by
    rw [R4C.stieltjes2C, R4C.trace_resolvC_sq hW hz.ne']
  rw [hval]
  refine hkey.congr_of_eventuallyEq ?_
  filter_upwards [hopen.mem_nhds hmem] with ζ hζ
  show R4C.stieltjesC W ζ = _
  rw [R4C.stieltjesC, R4C.trace_resolvC hW (ne_of_gt hζ)]

/-- `quad c z w = z (w - m)(w - m̃)` with `m̃ = (z m)⁻¹` the second root, for any root `m`.

Copy of `R1.quad_factor` (`RMT/R1.lean:347`); choice 8 keeps `RMT/General/` free of that
file, item F34. -/
private theorem quad_factor {z m : ℂ} (hz : z ≠ 0) (hmroot : MP.quad c z m = 0) (w : ℂ) :
    MP.quad c z w = z * (w - m) * (w - (z * m)⁻¹) := by
  have hm0 : m ≠ 0 := MP.root_ne_zero hmroot
  have hzm : z * m ≠ 0 := mul_ne_zero hz hm0
  have h : z * m ^ 2 + (z + 1 - (c : ℂ)) * m + 1 = 0 := hmroot
  set e : ℂ := (z * m)⁻¹ with he_def
  have he : z * m * e = 1 := mul_inv_cancel₀ hzm
  have hbr : z * m + z * e + z + 1 - (c : ℂ) = 0 := by
    have hmul : m * (z * m + z * e + z + 1 - (c : ℂ)) = 0 := by linear_combination h + he
    exact (mul_eq_zero.mp hmul).resolve_left hm0
  change z * w ^ 2 + (z + 1 - (c : ℂ)) * w + 1 = z * (w - m) * (w - e)
  linear_combination w * hbr - he

/-- Copy of `R1.norm_ge_rootLb` (`RMT/R1.lean:360`); choice 8 keeps `RMT/General/` free of
that file, item F34. -/
theorem norm_ge_rootLb (hz : z ≠ 0) {x r : ℂ}
    (hr : MP.quad c z x = -r) (hsmall : ‖r‖ ≤ 1 / 2) : rootLb c z ≤ ‖x‖ := by
  by_contra hcon
  rw [not_le] at hcon
  have hzn : 0 < ‖z‖ := norm_pos_iff.mpr hz
  have hD : 0 < ‖z‖ + ‖z + 1 - (c : ℂ)‖ := by positivity
  have hx1 : ‖x‖ ≤ 1 := le_of_lt (lt_of_lt_of_le hcon (min_le_left _ _))
  have hx2 : ‖x‖ < 1 / (2 * (‖z‖ + ‖z + 1 - (c : ℂ)‖)) :=
    lt_of_lt_of_le hcon (min_le_right _ _)
  have hxnn : (0 : ℝ) ≤ ‖x‖ := norm_nonneg _
  have h1 : -r - z * x ^ 2 - (z + 1 - (c : ℂ)) * x = 1 := by
    have h : z * x ^ 2 + (z + 1 - (c : ℂ)) * x + 1 = -r := hr
    linear_combination -h
  have hnorm : (1 : ℝ) ≤ ‖r‖ + ‖z‖ * ‖x‖ ^ 2 + ‖z + 1 - (c : ℂ)‖ * ‖x‖ := by
    have hone : (1 : ℝ) = ‖-r - z * x ^ 2 - (z + 1 - (c : ℂ)) * x‖ := by
      rw [h1]; simp
    have hstep : ‖-r - z * x ^ 2 - (z + 1 - (c : ℂ)) * x‖
        ≤ ‖-r‖ + ‖z * x ^ 2‖ + ‖(z + 1 - (c : ℂ)) * x‖ := by
      have e1 := norm_sub_le (-r - z * x ^ 2) ((z + 1 - (c : ℂ)) * x)
      have e2 := norm_sub_le (-r) (z * x ^ 2)
      linarith
    rw [hone]
    simpa [norm_mul, norm_pow] using hstep
  have hsqle : ‖x‖ ^ 2 ≤ ‖x‖ := by nlinarith
  have hlin : ‖z‖ * ‖x‖ ^ 2 + ‖z + 1 - (c : ℂ)‖ * ‖x‖
      ≤ (‖z‖ + ‖z + 1 - (c : ℂ)‖) * ‖x‖ := by
    nlinarith [mul_le_mul_of_nonneg_left hsqle hzn.le]
  have hhalf : (‖z‖ + ‖z + 1 - (c : ℂ)‖) * ‖x‖ < 1 / 2 := by
    have := (mul_lt_mul_of_pos_left hx2 hD)
    rw [mul_one_div] at this
    calc (‖z‖ + ‖z + 1 - (c : ℂ)‖) * ‖x‖
        < (‖z‖ + ‖z + 1 - (c : ℂ)‖) / (2 * (‖z‖ + ‖z + 1 - (c : ℂ)‖)) := this
      _ = 1 / 2 := by field_simp
  linarith

/-- Copy of `R1.norm_sub_root_le` (`RMT/R1.lean:395`); choice 8 keeps `RMT/General/` free of
that file, item F34. -/
theorem norm_sub_root_le (hc : 0 < c) (hz : 0 < z.im) {x r m : ℂ}
    (hm : 0 < m.im) (hmroot : MP.quad c z m = 0) (hx : 0 < x.im)
    (hr : MP.quad c z x = -r) : ‖x - m‖ ≤ ‖r‖ / (‖z‖ * x.im) := by
  have hz0 : z ≠ 0 := MP.ne_zero_of_im_pos hz
  have hzn : 0 < ‖z‖ := norm_pos_iff.mpr hz0
  have hm0 : m ≠ 0 := MP.root_ne_zero hmroot
  have hzm : z * m ≠ 0 := mul_ne_zero hz0 hm0
  have h : z * m ^ 2 + (z + 1 - (c : ℂ)) * m + 1 = 0 := hmroot
  set mt : ℂ := (z * m)⁻¹ with hmtdef
  have hmtne : mt ≠ 0 := inv_ne_zero hzm
  have hminv : m * m⁻¹ = 1 := mul_inv_cancel₀ hm0
  have hPQ : (1 + m⁻¹) * (1 + mt⁻¹) = (c : ℂ) := by
    have hmtinv : mt⁻¹ = z * m := by rw [hmtdef, inv_inv]
    rw [hmtinv]
    refine mul_left_cancel₀ hm0 ?_
    linear_combination h + (1 + z * m) * hminv
  have hmtim : mt.im ≤ 0 := by
    by_contra hcon
    rw [not_le] at hcon
    obtain ⟨h1, -⟩ := MP.im_eq_zero_of_im_nonneg hc hm0 hmtne hPQ hm.le hcon.le
    linarith
  have himx : x.im ≤ ‖x - mt‖ := by
    have h1 : (x - mt).im = x.im - mt.im := by simp
    calc x.im ≤ (x - mt).im := by rw [h1]; linarith
      _ ≤ |(x - mt).im| := le_abs_self _
      _ ≤ ‖x - mt‖ := Complex.abs_im_le_norm _
  have hfac : MP.quad c z x = z * (x - m) * (x - mt) := quad_factor hz0 hmroot x
  have hnr : ‖r‖ = ‖z‖ * ‖x - m‖ * ‖x - mt‖ := by
    have h0 : ‖r‖ = ‖MP.quad c z x‖ := by rw [hr, norm_neg]
    rw [h0, hfac, norm_mul, norm_mul]
  rw [le_div_iff₀ (by positivity : (0 : ℝ) < ‖z‖ * x.im), hnr]
  have hxm : (0 : ℝ) ≤ ‖x - m‖ := norm_nonneg _
  have key : ‖x - m‖ * x.im ≤ ‖x - m‖ * ‖x - mt‖ := mul_le_mul_of_nonneg_left himx hxm
  calc ‖x - m‖ * (‖z‖ * x.im) = ‖z‖ * (‖x - m‖ * x.im) := by ring
    _ ≤ ‖z‖ * (‖x - m‖ * ‖x - mt‖) := mul_le_mul_of_nonneg_left key hzn.le
    _ = ‖z‖ * ‖x - m‖ * ‖x - mt‖ := by ring

/-- A point of the closed disc of radius `ρ ≤ η` around `z` keeps `Im ≥ η - ρ`.

Copy of `R1.im_ge_of_mem_closedBall` (`RMT/R1.lean:570`); choice 8 keeps `RMT/General/` free
of that file, item F34. -/
theorem im_ge_of_mem_closedBall {ζ : ℂ} {ρ : ℝ} (hζ : ζ ∈ Metric.closedBall z ρ) :
    z.im - ρ ≤ ζ.im := by
  have h1 : ‖ζ - z‖ ≤ ρ := by
    rw [Metric.mem_closedBall, Complex.dist_eq] at hζ
    exact hζ
  have h2 : |(ζ - z).im| ≤ ‖ζ - z‖ := Complex.abs_im_le_norm _
  have h3 : (ζ - z).im = ζ.im - z.im := by simp
  rw [h3] at h2
  have h4 := abs_le.mp (h2.trans h1)
  linarith [h4.1]

/-- Copy of `R1.hasDerivAt_diff` (`RMT/R1.lean:581`); choice 8 keeps `RMT/General/` free of
that file, item F34. -/
theorem hasDerivAt_diff (hc : 0 < c) (_hd : 0 < d) (hW : W.IsHermitian) {ζ : ℂ}
    (hζ : 0 < ζ.im) :
    HasDerivAt (fun w => R4C.stieltjesC W w - MP.mC c w)
      (R4C.stieltjes2C W ζ - MP.mCDeriv c ζ) ζ :=
  (hasDerivAt_stieltjesC hW hζ).sub (MP.hasDerivAt_mC hc.le hζ)

/-- Copy of `R1.norm_deriv_diff_le` (`RMT/R1.lean:587`); choice 8 keeps `RMT/General/` free
of that file, item F34. -/
theorem norm_deriv_diff_le (hc : 0 < c) (hd : 0 < d) (hW : W.IsHermitian) {ζ : ℂ}
    (hζ : 0 < ζ.im) :
    ‖R4C.stieltjes2C W ζ - MP.mCDeriv c ζ‖ ≤ 2 / ζ.im ^ 2 := by
  have h1 := R4C.norm_stieltjes2C_le hW hζ hd
  have h2 := MP.norm_mCDeriv_le' hc hζ
  calc ‖R4C.stieltjes2C W ζ - MP.mCDeriv c ζ‖
      ≤ ‖R4C.stieltjes2C W ζ‖ + ‖MP.mCDeriv c ζ‖ := norm_sub_le _ _
    _ ≤ 1 / ζ.im ^ 2 + 1 / ζ.im ^ 2 := add_le_add h1 h2
    _ = 2 / ζ.im ^ 2 := by ring

/-- **Cauchy's estimate at radius `η/2`.** A uniform bound `C` on the circle bounds the gap of
the derivatives at the center by `C / (η/2)`.

Copy of `R1.norm_stieltjes2C_sub_mCDeriv_le` (`RMT/R1.lean:599`); choice 8 keeps
`RMT/General/` free of that file, item F34. -/
theorem norm_stieltjes2C_sub_mCDeriv_le (hc : 0 < c) (hz : 0 < z.im) (hd : 0 < d)
    (hW : W.IsHermitian) {C : ℝ}
    (hC : ∀ ζ ∈ Metric.sphere z (z.im / 2), ‖R4C.stieltjesC W ζ - MP.mC c ζ‖ ≤ C) :
    ‖R4C.stieltjes2C W z - MP.mCDeriv c z‖ ≤ C / (z.im / 2) := by
  have hrpos : (0 : ℝ) < z.im / 2 := by linarith
  have hball : ∀ ζ ∈ Metric.closedBall z (z.im / 2), 0 < ζ.im := by
    intro ζ hζ
    have := im_ge_of_mem_closedBall (z := z) hζ
    linarith
  have hdiffOn : DifferentiableOn ℂ (fun w => R4C.stieltjesC W w - MP.mC c w)
      (Metric.closedBall z (z.im / 2)) := fun ζ hζ =>
    ((hasDerivAt_diff hc hd hW (hball ζ hζ)).differentiableAt).differentiableWithinAt
  have hdc : DiffContOnCl ℂ (fun w => R4C.stieltjesC W w - MP.mC c w)
      (Metric.ball z (z.im / 2)) := by
    constructor
    · exact hdiffOn.mono Metric.ball_subset_closedBall
    · rw [closure_ball z (ne_of_gt hrpos)]
      exact hdiffOn.continuousOn
  have hderiv : deriv (fun w => R4C.stieltjesC W w - MP.mC c w) z
      = R4C.stieltjes2C W z - MP.mCDeriv c z := (hasDerivAt_diff hc hd hW hz).deriv
  rw [← hderiv]
  exact Complex.norm_deriv_le_of_forall_mem_sphere_norm_le hrpos hdc hC

/-- **Equicontinuity in `ζ`.** On the closed disc of radius `3η/4` the gap `s_N - m` is
`32/η²`-Lipschitz, uniformly in `W`.

Copy of `R1.norm_diff_sub_diff_le` (`RMT/R1.lean:624`); choice 8 keeps `RMT/General/` free of
that file, item F34. -/
theorem norm_diff_sub_diff_le (hc : 0 < c) (hz : 0 < z.im) (hd : 0 < d) (hW : W.IsHermitian)
    {ζ ζ' : ℂ} (hζ : ζ ∈ Metric.closedBall z (3 * z.im / 4))
    (hζ' : ζ' ∈ Metric.closedBall z (3 * z.im / 4)) :
    ‖(R4C.stieltjesC W ζ - MP.mC c ζ) - (R4C.stieltjesC W ζ' - MP.mC c ζ')‖
      ≤ (32 / z.im ^ 2) * ‖ζ - ζ'‖ := by
  have him : ∀ w ∈ Metric.closedBall z (3 * z.im / 4), z.im / 4 ≤ w.im := by
    intro w hw
    have := im_ge_of_mem_closedBall (z := z) hw
    linarith
  have hbound : ∀ w ∈ Metric.closedBall z (3 * z.im / 4),
      ‖R4C.stieltjes2C W w - MP.mCDeriv c w‖ ≤ 32 / z.im ^ 2 := by
    intro w hw
    have h1 : z.im / 4 ≤ w.im := him w hw
    have hwpos : 0 < w.im := by linarith
    refine (norm_deriv_diff_le hc hd hW hwpos).trans ?_
    rw [div_le_div_iff₀ (by positivity) (by positivity)]
    nlinarith [sq_nonneg (w.im - z.im / 4), sq_nonneg z.im]
  have hder : ∀ w ∈ Metric.closedBall z (3 * z.im / 4),
      HasDerivWithinAt (fun v => R4C.stieltjesC W v - MP.mC c v)
        (R4C.stieltjes2C W w - MP.mCDeriv c w) (Metric.closedBall z (3 * z.im / 4)) w := by
    intro w hw
    have h1 : z.im / 4 ≤ w.im := him w hw
    exact (hasDerivAt_diff hc hd hW (by linarith)).hasDerivWithinAt
  exact (convex_closedBall z (3 * z.im / 4)).norm_image_sub_le_of_norm_hasDerivWithin_le
    hder hbound hζ' hζ

end GenRMT
end StackedSVD
