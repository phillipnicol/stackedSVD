/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.R0
import StackedSVD.RMT.R4C
import StackedSVD.RMT.ResolvDeriv
import StackedSVD.RMT.SteinStep
import StackedSVD.RMT.MP7
import StackedSVD.Prob.GaussianAdapters

/-!
# Item R1: the complex Marchenko-Pastur law for the trace of the resolvent

Specification: `notes/archive/rmt_R1.md` (restructured 2026-08-30, decision D7, choices 27 to 31).

Setting. `Y ~ gaussianMatrix p d` is the noise block `B` of item R0, `W₀ = d⁻¹ Yᵀ Y`,
`G(z) = (W₀ - z)⁻¹`, `s_N(z) = d⁻¹ tr G(z)` and `c_N = p / d`. For `Im z > 0` the two
statements are

* **R1a** `s_N(z) → MP.mC c z` in probability,
* **R1b** `d⁻¹ tr G(z)² → MP.mCDeriv c z` in probability,

both exported in the norm form `TendstoInProb μ (fun N ω => ‖· - ·‖) 0` (choice 28), which
items R2, T and R3⁻ consume.

## Contents

1. `variance_le_of_lipschitz`: `Var F ≤ 4 L²` for a Lipschitz `F` under `gaussianMatrix`, by
   layer cake from the proved tail bound `measure_abs_ge_le_of_lipschitz`. The note asks for
   this lemma in `Prob/GaussianAdapters.lean`; it lives here so that R1 owns one file.
2. `im_stieltjesC_ge`, `hasDerivAt_stieltjesC`: the two deterministic facts of step 1 that
   `RMT/R4C.lean` leaves open.
3. `lipschitzWith_stieltjesC_re`, `lipschitzWith_stieltjesC_im`: the global Lipschitz bound of
   step 3, in the Frobenius metric through `matrixEquivE`.
4. `norm_quad_integral_le`: the Stein equation of step 2, in the residual form.
5. `norm_ge_rootLb`, `norm_sub_root_le`: the stability of the MP quadratic, step 4.
6. `tendstoInProb_stieltjesC` (R1a) and `tendstoInProb_stieltjes2C` (R1b).

No `axiom`, no `sorry`. Step 3 (`lipschitzWith_stieltjesC_re`, `_im`) closed 2026-08-30 from
`RMT/ResolvDeriv.lean`; step 2 (`norm_quad_integral_le`) closed 2026-08-30 from
`RMT/SteinStep.lean` (`notes/archive/agent_reports/proof_r1_stein.md`).
-/

open Filter Topology MeasureTheory ProbabilityTheory Set
open scoped Matrix ENNReal NNReal

namespace StackedSVD

/-! ### The variance bound from the Lipschitz tail (choice 31, decision D7) -/

/-- **`Var F ≤ 4 L²`.** For `F` `L`-Lipschitz on the flattened matrix space, the variance of
`F` under `gaussianMatrix p d` is at most `4 L²`. Proof: the layer cake
`Integrable.integral_eq_integral_meas_le` applied to `(F - E F)²`, whose tail at level `t` is
the tail of `|F - E F|` at level `√t`, then `measure_abs_ge_le_of_lipschitz` and
`∫₀^∞ 2 exp(-t/(2L²)) dt = 4 L²`.

`notes/archive/rmt_R1.md` asks for this lemma in `Prob/GaussianAdapters.lean`; the name is the same.
-/
theorem variance_le_of_lipschitz {p d : ℕ} {LL : ℝ} (hpd : 0 < p * d) (hL : 0 < LL)
    {F : EuclideanSpace ℝ (Fin (p * d)) → ℝ} (hF : LipschitzWith ⟨LL, hL.le⟩ F) :
    ∫ Z, (F (matrixEquivE p d Z) - ∫ Z', F (matrixEquivE p d Z') ∂gaussianMatrix p d) ^ 2
        ∂gaussianMatrix p d
      ≤ 4 * LL ^ 2 := by
  have hLnn : (0 : ℝ≥0) < (⟨LL, hL.le⟩ : ℝ≥0) := NNReal.coe_pos.mp hL
  set cst : ℝ := ∫ Z', F (matrixEquivE p d Z') ∂gaussianMatrix p d with hcst
  set g : Matrix (Fin p) (Fin d) ℝ → ℝ :=
    fun Z => (F (matrixEquivE p d Z) - cst) ^ 2 with hgdef
  have hmeasF : Measurable fun Z => F (matrixEquivE p d Z) :=
    hF.continuous.measurable.comp (matrixEquivE p d).measurable
  have hgmeas : Measurable g := (hmeasF.sub measurable_const).pow_const 2
  have hgnn : 0 ≤ᵐ[gaussianMatrix p d] g := Filter.Eventually.of_forall fun _ => sq_nonneg _
  set a : ℝ := -(2 * LL ^ 2)⁻¹ with hadef
  have ha : a < 0 := by
    rw [hadef, neg_lt_zero]
    positivity
  -- the tail of `g` at level `t` is the tail of `|F - E F|` at level `√t`
  have hbound : ∀ t : ℝ, 0 < t →
      gaussianMatrix p d {Z | t ≤ g Z} ≤ ENNReal.ofReal (2 * Real.exp (a * t)) := by
    intro t ht
    have hst : 0 < Real.sqrt t := Real.sqrt_pos.mpr ht
    have hset : {Z | t ≤ g Z}
        = {Z | Real.sqrt t ≤ |F (matrixEquivE p d Z) - cst|} := by
      ext Z
      simp only [Set.mem_ofPred_eq, hgdef]
      constructor
      · intro h
        have h1 : Real.sqrt t ≤ Real.sqrt ((F (matrixEquivE p d Z) - cst) ^ 2) :=
          Real.sqrt_le_sqrt h
        rwa [Real.sqrt_sq_eq_abs] at h1
      · intro h
        have h2 : Real.sqrt t ^ 2 ≤ |F (matrixEquivE p d Z) - cst| ^ 2 :=
          pow_le_pow_left₀ (Real.sqrt_nonneg t) h 2
        rwa [Real.sq_sqrt ht.le, sq_abs] at h2
    have hkey := measure_abs_ge_le_of_lipschitz (p := p) (d := d) hpd hLnn hF hst
    rw [← hcst] at hkey
    have hsq : Real.sqrt t ^ 2 = t := Real.sq_sqrt ht.le
    rw [hsq] at hkey
    have hkey2 :
        ((gaussianMatrix p d) {Z | Real.sqrt t ≤ |F (matrixEquivE p d Z) - cst|}).toReal
          ≤ 2 * Real.exp (-t / (2 * LL ^ 2)) := hkey
    have harg : -t / (2 * LL ^ 2) = a * t := by
      rw [hadef]; ring
    rw [harg] at hkey2
    have hfin : gaussianMatrix p d {Z | Real.sqrt t ≤ |F (matrixEquivE p d Z) - cst|} ≠ ⊤ :=
      measure_ne_top _ _
    rw [hset, ← ENNReal.ofReal_toReal hfin]
    exact ENNReal.ofReal_le_ofReal hkey2
  -- layer cake
  have hlc : ∫⁻ Z, ENNReal.ofReal (g Z) ∂(gaussianMatrix p d)
      = ∫⁻ t in Ioi 0, gaussianMatrix p d {Z | t ≤ g Z} :=
    lintegral_eq_lintegral_meas_le _ hgnn hgmeas.aemeasurable
  have hmono : ∫⁻ t in Ioi 0, gaussianMatrix p d {Z | t ≤ g Z}
      ≤ ∫⁻ t in Ioi 0, ENNReal.ofReal (2 * Real.exp (a * t)) := by
    refine lintegral_mono_ae ?_
    filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
    exact hbound t ht
  have hexp : ∫⁻ t in Ioi 0, ENNReal.ofReal (2 * Real.exp (a * t))
      = ENNReal.ofReal (4 * LL ^ 2) := by
    have hint : IntegrableOn (fun t : ℝ => 2 * Real.exp (a * t)) (Ioi 0) :=
      (integrableOn_exp_mul_Ioi ha 0).const_mul 2
    have hnn : 0 ≤ᵐ[volume.restrict (Ioi 0)] fun t : ℝ => 2 * Real.exp (a * t) :=
      Filter.Eventually.of_forall fun t => by positivity
    rw [← ofReal_integral_eq_lintegral_ofReal hint hnn]
    congr 1
    rw [integral_const_mul, integral_exp_mul_Ioi ha 0]
    rw [hadef]
    simp only [mul_zero, Real.exp_zero]
    field_simp
    norm_num
  have hfinal : ∫⁻ Z, ENNReal.ofReal (g Z) ∂(gaussianMatrix p d)
      ≤ ENNReal.ofReal (4 * LL ^ 2) := by
    rw [hlc]
    exact hmono.trans (le_of_eq hexp)
  have hval : ∫ Z, g Z ∂(gaussianMatrix p d)
      = (∫⁻ Z, ENNReal.ofReal (g Z) ∂(gaussianMatrix p d)).toReal :=
    integral_eq_lintegral_of_nonneg_ae hgnn hgmeas.aestronglyMeasurable
  change ∫ Z, g Z ∂(gaussianMatrix p d) ≤ 4 * LL ^ 2
  rw [hval]
  calc (∫⁻ Z, ENNReal.ofReal (g Z) ∂(gaussianMatrix p d)).toReal
      ≤ (ENNReal.ofReal (4 * LL ^ 2)).toReal :=
        ENNReal.toReal_mono ENNReal.ofReal_ne_top hfinal
    _ = 4 * LL ^ 2 := ENNReal.toReal_ofReal (by positivity)

namespace R1

variable {p d : ℕ}

/-- The Wishart block of `Y = √d E⊥`. The complex resolvent and the two traces are proved in
`RMT/R4C.lean`. -/
noncomputable def W0 (Y : Matrix (Fin p) (Fin d) ℝ) : Matrix (Fin d) (Fin d) ℝ :=
  (d : ℝ)⁻¹ • (Yᵀ * Y)

/-- The Lipschitz constant of step 3. -/
noncomputable def lipConst (z : ℂ) (d : ℕ) : ℝ := 2 * Real.sqrt (z.im + ‖z‖) / (d * z.im ^ 2)

/-- The residual constant of step 2. -/
noncomputable def resConst (z : ℂ) : ℝ :=
  1 / z.im + ‖z‖ / z.im ^ 2 + 32 * ‖z‖ * (z.im + ‖z‖) / z.im ^ 4

/-- The `δ` of step 4, a lower bound on `‖E s_N‖`. -/
noncomputable def rootLb (c : ℝ) (z : ℂ) : ℝ := min 1 (1 / (2 * (‖z‖ + ‖z + 1 - c‖)))

theorem lipConst_nonneg (z : ℂ) (d : ℕ) : 0 ≤ lipConst z d := by
  have h : 0 ≤ z.im + ‖z‖ := by
    have := Complex.abs_im_le_norm z
    rcases abs_cases z.im with h1 | h1 <;> linarith [h1.1, h1.2]
  unfold lipConst
  positivity

/-- `lipConst` as an `ℝ≥0`, the shape `LipschitzWith` takes. The note writes
`⟨lipConst z d, sorry⟩`; this is the same term with the nonnegativity proof filled in. -/
noncomputable def lipConstNN (z : ℂ) (d : ℕ) : ℝ≥0 := ⟨lipConst z d, lipConst_nonneg z d⟩

@[simp]
theorem coe_lipConstNN (z : ℂ) (d : ℕ) : ((lipConstNN z d : ℝ≥0) : ℝ) = lipConst z d := rfl

variable {W : Matrix (Fin d) (Fin d) ℝ} {z : ℂ} {c : ℝ}

/-! ### Step 1: the two deterministic facts that `RMT/R4C.lean` leaves open -/

/-- `Im s ≥ η ‖s‖²`, by Cauchy-Schwarz on the eigenvalue sum. -/
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

/-- `d⁻¹ tr G(z)²` is the derivative in `z` of `d⁻¹ tr G(z)`. -/
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

/-! ### Steps 2 and 3: Stein and the Lipschitz bound -/

/-- Step 3, real part. `∇_Y Re s_N = -(2/d²) Re (Y G²)` has Frobenius norm at most
`lipConst z d`. -/
theorem lipschitzWith_stieltjesC_re (hz : 0 < z.im) (hd : 0 < d) :
    LipschitzWith (lipConstNN z d)
      (fun x : EuclideanSpace ℝ (Fin (p * d)) =>
        (R4C.stieltjesC (W0 ((matrixEquivE p d).symm x)) z).re) := by
  have hfun : (fun x : EuclideanSpace ℝ (Fin (p * d)) =>
        (R4C.stieltjesC (W0 ((matrixEquivE p d).symm x)) z).re)
      = fun x : EuclideanSpace ℝ (Fin (p * d)) =>
        (ResolvDeriv.sfun (p := p) (d := d) z x).re := by
    funext x
    rw [ResolvDeriv.sfun_eq]
    rfl
  rw [hfun]
  exact ResolvDeriv.lipschitzWith_sfun_re hz hd le_rfl

/-- Step 3, imaginary part. -/
theorem lipschitzWith_stieltjesC_im (hz : 0 < z.im) (hd : 0 < d) :
    LipschitzWith (lipConstNN z d)
      (fun x : EuclideanSpace ℝ (Fin (p * d)) =>
        (R4C.stieltjesC (W0 ((matrixEquivE p d).symm x)) z).im) := by
  have hfun : (fun x : EuclideanSpace ℝ (Fin (p * d)) =>
        (R4C.stieltjesC (W0 ((matrixEquivE p d).symm x)) z).im)
      = fun x : EuclideanSpace ℝ (Fin (p * d)) =>
        (ResolvDeriv.sfun (p := p) (d := d) z x).im := by
    funext x
    rw [ResolvDeriv.sfun_eq]
    rfl
  rw [hfun]
  exact ResolvDeriv.lipschitzWith_sfun_im hz hd le_rfl

/-- Step 2. The mean nearly solves the MP quadratic at `c_N = p / d`.

`RMT/SteinStep.lean` proves the Gaussian integration by parts and the residual form
`SteinStep.norm_quad_integral_le_aux`, which takes the two variance bounds of choice 31 as
hypotheses. This proof supplies them from `variance_le_of_lipschitz` at `L = lipConst z d`
and then compares `32 ‖z‖ (η + ‖z‖) / (d² η⁴)` with the third term of `resConst z / d`. -/
theorem norm_quad_integral_le (hz : 0 < z.im) (hp : 0 < p) (hd : 0 < d) :
    ‖MP.quad ((p : ℝ) / d) z (∫ Y, R4C.stieltjesC (W0 Y) z ∂gaussianMatrix p d)‖
      ≤ resConst z / d := by
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd
  have hd1 : (1 : ℝ) ≤ (d : ℝ) := by exact_mod_cast hd
  have hpd : 0 < p * d := Nat.mul_pos hp hd
  have hLpos : 0 < lipConst z d := by
    have h1 : 0 < z.im + ‖z‖ := by linarith [norm_nonneg z]
    have h2 : 0 < Real.sqrt (z.im + ‖z‖) := Real.sqrt_pos.mpr h1
    unfold lipConst
    positivity
  have hnn : 0 ≤ z.im + ‖z‖ := by positivity
  have hznn : (0 : ℝ) ≤ ‖z‖ := norm_nonneg z
  -- the two variance bounds of choice 31
  have hre : ∫ Z, ((ResolvDeriv.sfun z (matrixEquivE p d Z)).re
      - ∫ Z', (ResolvDeriv.sfun z (matrixEquivE p d Z')).re ∂gaussianMatrix p d) ^ 2
      ∂gaussianMatrix p d ≤ 4 * lipConst z d ^ 2 :=
    variance_le_of_lipschitz hpd hLpos
      (ResolvDeriv.lipschitzWith_sfun_re (p := p) (d := d) hz hd
        (C := ⟨lipConst z d, hLpos.le⟩) le_rfl)
  have him : ∫ Z, ((ResolvDeriv.sfun z (matrixEquivE p d Z)).im
      - ∫ Z', (ResolvDeriv.sfun z (matrixEquivE p d Z')).im ∂gaussianMatrix p d) ^ 2
      ∂gaussianMatrix p d ≤ 4 * lipConst z d ^ 2 :=
    variance_le_of_lipschitz hpd hLpos
      (ResolvDeriv.lipschitzWith_sfun_im (p := p) (d := d) hz hd
        (C := ⟨lipConst z d, hLpos.le⟩) le_rfl)
  have hkey := SteinStep.norm_quad_integral_le_aux (p := p) (d := d) (z := z) hz hd hre him
  refine hkey.trans ?_
  -- the arithmetic of the three residual terms
  have hsq : lipConst z d ^ 2 = 4 * (z.im + ‖z‖) / ((d : ℝ) ^ 2 * z.im ^ 4) := by
    unfold lipConst
    rw [div_pow, mul_pow, Real.sq_sqrt hnn, mul_pow]
    ring
  have hLHS : ‖z‖ * (2 * (4 * lipConst z d ^ 2))
      = 32 * ‖z‖ * (z.im + ‖z‖) / ((d : ℝ) ^ 2 * z.im ^ 4) := by
    rw [hsq]
    field_simp
    ring
  have hcmp : 32 * ‖z‖ * (z.im + ‖z‖) / ((d : ℝ) ^ 2 * z.im ^ 4)
      ≤ 32 * ‖z‖ * (z.im + ‖z‖) / ((d : ℝ) * z.im ^ 4) := by
    gcongr
    all_goals first
      | positivity
      | nlinarith [pow_pos hz 4]
  have hres : resConst z / (d : ℝ)
      = 1 / ((d : ℝ) * z.im) + ‖z‖ / ((d : ℝ) * z.im ^ 2)
        + 32 * ‖z‖ * (z.im + ‖z‖) / ((d : ℝ) * z.im ^ 4) := by
    unfold resConst
    field_simp
  have hgoal : resConst z / ((d : ℕ) : ℝ) = resConst z / (d : ℝ) := rfl
  rw [hgoal, hres, hLHS]
  linarith

/-! ### Step 4: root stability, pure algebra on `MP.quad` -/

/-- `quad c z w = z (w - m)(w - m̃)` with `m̃ = (z m)⁻¹` the second root, for any root `m`. -/
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

/-! ### Bounds, measurability and integrability of `Y ↦ s_N(z)`

Every statement here is deterministic in `Y` or a routine consequence of the bound
`‖s_N(z)‖ ≤ 1/η` of `RMT/R4C.lean`. -/

theorem isHermitian_W0 (Y : Matrix (Fin p) (Fin d) ℝ) : (W0 Y).IsHermitian := by
  have h := isHermitian_transpose_mul_self Y
  change ((d : ℝ)⁻¹ • (Yᵀ * Y))ᴴ = (d : ℝ)⁻¹ • (Yᵀ * Y)
  rw [Matrix.conjTranspose_smul, star_trivial, h.eq]

theorem norm_stieltjesC_W0_le (hz : 0 < z.im) (hd : 0 < d) (Y : Matrix (Fin p) (Fin d) ℝ) :
    ‖R4C.stieltjesC (W0 Y) z‖ ≤ 1 / z.im :=
  R4C.norm_stieltjesC_le (isHermitian_W0 Y) hz hd

theorem rootLb_pos (hz : z ≠ 0) : 0 < rootLb c z := by
  have hzn : 0 < ‖z‖ := norm_pos_iff.mpr hz
  have hnn : (0 : ℝ) ≤ ‖z + 1 - (c : ℂ)‖ := norm_nonneg _
  refine lt_min one_pos ?_
  positivity

theorem lipConst_pos (hz : 0 < z.im) (hd : 0 < d) : 0 < lipConst z d := by
  have h1 : 0 < z.im + ‖z‖ := by linarith [norm_nonneg z]
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd
  have h2 : 0 < Real.sqrt (z.im + ‖z‖) := Real.sqrt_pos.mpr h1
  unfold lipConst
  positivity

theorem measurable_stieltjesC_re (hz : 0 < z.im) (hd : 0 < d) :
    Measurable fun Y : Matrix (Fin p) (Fin d) ℝ => (R4C.stieltjesC (W0 Y) z).re := by
  have h := ((lipschitzWith_stieltjesC_re (p := p) (d := d) hz hd).continuous.measurable).comp
    (matrixEquivE p d).measurable
  simpa [Function.comp_def] using h

theorem measurable_stieltjesC_im (hz : 0 < z.im) (hd : 0 < d) :
    Measurable fun Y : Matrix (Fin p) (Fin d) ℝ => (R4C.stieltjesC (W0 Y) z).im := by
  have h := ((lipschitzWith_stieltjesC_im (p := p) (d := d) hz hd).continuous.measurable).comp
    (matrixEquivE p d).measurable
  simpa [Function.comp_def] using h

theorem measurable_stieltjesC (hz : 0 < z.im) (hd : 0 < d) :
    Measurable fun Y : Matrix (Fin p) (Fin d) ℝ => R4C.stieltjesC (W0 Y) z := by
  have hre := measurable_stieltjesC_re (p := p) (d := d) hz hd
  have him := measurable_stieltjesC_im (p := p) (d := d) hz hd
  have hsplit : (fun Y : Matrix (Fin p) (Fin d) ℝ => R4C.stieltjesC (W0 Y) z)
      = fun Y => (((R4C.stieltjesC (W0 Y) z).re : ℝ) : ℂ)
        + (((R4C.stieltjesC (W0 Y) z).im : ℝ) : ℂ) * Complex.I := by
    funext Y
    exact (Complex.re_add_im _).symm
  rw [hsplit]
  exact (Complex.measurable_ofReal.comp hre).add
    ((Complex.measurable_ofReal.comp him).mul measurable_const)

theorem integrable_stieltjesC (hz : 0 < z.im) (hd : 0 < d) :
    Integrable (fun Y : Matrix (Fin p) (Fin d) ℝ => R4C.stieltjesC (W0 Y) z)
      (gaussianMatrix p d) :=
  Integrable.mono' (integrable_const (1 / z.im))
    (measurable_stieltjesC hz hd).aestronglyMeasurable
    (Filter.Eventually.of_forall fun Y => norm_stieltjesC_W0_le hz hd Y)

/-- Jensen for a bounded measurable real function under a probability measure. -/
private theorem sq_integral_le {α : Type*} [MeasurableSpace α] (ν : Measure α)
    [IsProbabilityMeasure ν] {f : α → ℝ} {C : ℝ} (hf : Measurable f) (hb : ∀ x, |f x| ≤ C) :
    (∫ x, f x ∂ν) ^ 2 ≤ ∫ x, f x ^ 2 ∂ν := by
  have hmem : MemLp f 2 ν :=
    MemLp.mono_exponent (memLp_top_of_bound hf.aestronglyMeasurable C
      (Filter.Eventually.of_forall fun x => by simpa [Real.norm_eq_abs] using hb x)) le_top
  have h := variance_nonneg f ν
  rw [variance_eq_sub hmem] at h
  simp only [Pi.pow_apply] at h
  linarith

/-- Step 1 under the integral: `Im E s ≥ η ‖E s‖²`, by `im_stieltjesC_ge` and Jensen. -/
theorem im_integral_stieltjesC_ge (hz : 0 < z.im) (hd : 0 < d) :
    z.im * ‖∫ Y, R4C.stieltjesC (W0 Y) z ∂gaussianMatrix p d‖ ^ 2
      ≤ (∫ Y, R4C.stieltjesC (W0 Y) z ∂gaussianMatrix p d).im := by
  set ν : Measure (Matrix (Fin p) (Fin d) ℝ) := gaussianMatrix p d with hν
  set f : Matrix (Fin p) (Fin d) ℝ → ℂ := fun Y => R4C.stieltjesC (W0 Y) z with hfdef
  have hbd : ∀ Y, ‖f Y‖ ≤ 1 / z.im := fun Y => norm_stieltjesC_W0_le hz hd Y
  have hbre : ∀ Y, |(f Y).re| ≤ 1 / z.im := fun Y =>
    le_trans (Complex.abs_re_le_norm (f Y)) (hbd Y)
  have hbim : ∀ Y, |(f Y).im| ≤ 1 / z.im := fun Y =>
    le_trans (Complex.abs_im_le_norm (f Y)) (hbd Y)
  have hmre : Measurable fun Y => (f Y).re := measurable_stieltjesC_re hz hd
  have hmim : Measurable fun Y => (f Y).im := measurable_stieltjesC_im hz hd
  have hIc : Integrable f ν := integrable_stieltjesC hz hd
  have hIre : Integrable (fun Y => (f Y).re) ν :=
    Integrable.mono' (integrable_const (1 / z.im)) hmre.aestronglyMeasurable
      (Filter.Eventually.of_forall fun Y => by simpa [Real.norm_eq_abs] using hbre Y)
  have hIim : Integrable (fun Y => (f Y).im) ν :=
    Integrable.mono' (integrable_const (1 / z.im)) hmim.aestronglyMeasurable
      (Filter.Eventually.of_forall fun Y => by simpa [Real.norm_eq_abs] using hbim Y)
  have hIre2 : Integrable (fun Y => (f Y).re ^ 2) ν :=
    Integrable.mono' (integrable_const ((1 / z.im) ^ 2))
      (hmre.pow_const 2).aestronglyMeasurable
      (Filter.Eventually.of_forall fun Y => by
        rw [Real.norm_eq_abs, abs_pow]
        exact pow_le_pow_left₀ (abs_nonneg _) (hbre Y) 2)
  have hIim2 : Integrable (fun Y => (f Y).im ^ 2) ν :=
    Integrable.mono' (integrable_const ((1 / z.im) ^ 2))
      (hmim.pow_const 2).aestronglyMeasurable
      (Filter.Eventually.of_forall fun Y => by
        rw [Real.norm_eq_abs, abs_pow]
        exact pow_le_pow_left₀ (abs_nonneg _) (hbim Y) 2)
  set A : ℂ := ∫ Y, f Y ∂ν with hA
  have hAre : A.re = ∫ Y, (f Y).re ∂ν := by
    have h := ContinuousLinearMap.integral_comp_comm Complex.reCLM hIc
    simpa [hA] using h.symm
  have hAim : A.im = ∫ Y, (f Y).im ∂ν := by
    have h := ContinuousLinearMap.integral_comp_comm Complex.imCLM hIc
    simpa [hA] using h.symm
  have hnormsq : ‖A‖ ^ 2 = A.re ^ 2 + A.im ^ 2 := by
    rw [Complex.sq_norm, Complex.normSq_apply]; ring
  have hjre : A.re ^ 2 ≤ ∫ Y, (f Y).re ^ 2 ∂ν := by
    rw [hAre]; exact sq_integral_le ν hmre hbre
  have hjim : A.im ^ 2 ≤ ∫ Y, (f Y).im ^ 2 ∂ν := by
    rw [hAim]; exact sq_integral_le ν hmim hbim
  have hpoint : ∀ Y, z.im * ((f Y).re ^ 2 + (f Y).im ^ 2) ≤ (f Y).im := by
    intro Y
    have h := im_stieltjesC_ge (isHermitian_W0 Y) hz
    rw [Complex.sq_norm, Complex.normSq_apply] at h
    calc z.im * ((f Y).re ^ 2 + (f Y).im ^ 2)
        = z.im * ((f Y).re * (f Y).re + (f Y).im * (f Y).im) := by ring
      _ ≤ (f Y).im := h
  have hint : ∫ Y, z.im * ((f Y).re ^ 2 + (f Y).im ^ 2) ∂ν ≤ ∫ Y, (f Y).im ∂ν :=
    integral_mono ((hIre2.add hIim2).const_mul _) hIim hpoint
  have heq : ∫ Y, z.im * ((f Y).re ^ 2 + (f Y).im ^ 2) ∂ν
      = z.im * ((∫ Y, (f Y).re ^ 2 ∂ν) + ∫ Y, (f Y).im ^ 2 ∂ν) := by
    rw [integral_const_mul, integral_add hIre2 hIim2]
  rw [heq] at hint
  calc z.im * ‖A‖ ^ 2 = z.im * (A.re ^ 2 + A.im ^ 2) := by rw [hnormsq]
    _ ≤ z.im * ((∫ Y, (f Y).re ^ 2 ∂ν) + ∫ Y, (f Y).im ^ 2 ∂ν) := by
        exact mul_le_mul_of_nonneg_left (by linarith) hz.le
    _ ≤ ∫ Y, (f Y).im ∂ν := hint
    _ = A.im := hAim.symm

/-! ### Step 6: the Cauchy estimate and the equicontinuity in `z` -/

/-- A point of the closed disc of radius `ρ ≤ η` around `z` keeps `Im ≥ η - ρ`. -/
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

theorem hasDerivAt_diff (hc : 0 < c) (_hd : 0 < d) (hW : W.IsHermitian) {ζ : ℂ}
    (hζ : 0 < ζ.im) :
    HasDerivAt (fun w => R4C.stieltjesC W w - MP.mC c w)
      (R4C.stieltjes2C W ζ - MP.mCDeriv c ζ) ζ :=
  (hasDerivAt_stieltjesC hW hζ).sub (MP.hasDerivAt_mC hc.le hζ)

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
the derivatives at the center by `C / (η/2)`. -/
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
`32/η²`-Lipschitz, uniformly in `W`. -/
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

/-- Changing `c` in the quadratic costs `|c - c'| ‖w‖`. -/
private theorem quad_sub_quad (c c' : ℝ) (z w : ℂ) :
    MP.quad c z w - MP.quad c' z w = ((c' : ℂ) - (c : ℂ)) * w := by
  change (z * w ^ 2 + (z + 1 - (c : ℂ)) * w + 1) - (z * w ^ 2 + (z + 1 - (c' : ℂ)) * w + 1) = _
  ring

/-- **Step 4 in the limit.** The mean `E s_N` converges to `MP.mC c z`: the Stein residual of
step 2 plus the change of `c_N` to `c` is `o(1)`, and step 4 turns that into a distance to the
root. -/
theorem tendsto_norm_integral_sub_mC {pN dN : ℕ → ℕ} (hcpos : 0 < c) (hzim : 0 < z.im)
    (hdtop : Tendsto dN atTop atTop) (hppos : ∀ N, 0 < pN N)
    (hratio : Tendsto (fun N => (pN N : ℝ) / dN N) atTop (𝓝 c)) :
    Tendsto (fun N => ‖(∫ Yv, R4C.stieltjesC (W0 Yv) z ∂gaussianMatrix (pN N) (dN N))
        - MP.mC c z‖) atTop (𝓝 0) := by
  have hz0 : z ≠ 0 := MP.ne_zero_of_im_pos hzim
  have hzn : 0 < ‖z‖ := norm_pos_iff.mpr hz0
  have hδ : 0 < rootLb c z := rootLb_pos hz0
  have hDen : 0 < ‖z‖ * (z.im * rootLb c z ^ 2) := by positivity
  set E : ℕ → ℂ :=
    fun N => ∫ Yv, R4C.stieltjesC (W0 Yv) z ∂gaussianMatrix (pN N) (dN N) with hEdef
  set R : ℕ → ℝ := fun N => ‖MP.quad c z (E N)‖ with hRdef
  have hev1 : ∀ᶠ N in atTop, 0 < dN N := hdtop.eventually_gt_atTop 0
  have hdRto : Tendsto (fun N => ((dN N : ℝ))) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp hdtop
  -- the bounding sequence of the residual
  set G : ℕ → ℝ :=
    fun N => resConst z / (dN N : ℝ) + |(pN N : ℝ) / dN N - c| * (1 / z.im) with hGdef
  have hGto : Tendsto G atTop (𝓝 0) := by
    have h1 : Tendsto (fun N => resConst z / (dN N : ℝ)) atTop (𝓝 0) :=
      Filter.Tendsto.div_atTop tendsto_const_nhds hdRto
    have h2 : Tendsto (fun N => |(pN N : ℝ) / dN N - c|) atTop (𝓝 0) := by
      have := (hratio.sub_const c).abs
      simpa using this
    have h3 : Tendsto (fun N => |(pN N : ℝ) / dN N - c| * (1 / z.im)) atTop (𝓝 0) := by
      simpa using h2.mul_const (1 / z.im)
    simpa [hGdef] using h1.add h3
  have hRle : ∀ᶠ N in atTop, R N ≤ G N := by
    filter_upwards [hev1] with N hdN
    have hpN := hppos N
    have h1 := norm_quad_integral_le (p := pN N) (d := dN N) hzim hpN hdN
    have hEnorm : ‖E N‖ ≤ 1 / z.im := by
      have h := norm_integral_le_of_norm_le_const
        (μ := gaussianMatrix (pN N) (dN N)) (C := 1 / z.im)
        (Filter.Eventually.of_forall fun Yv => norm_stieltjesC_W0_le hzim hdN Yv)
      simpa [hEdef] using h
    have hsub : MP.quad c z (E N) - MP.quad ((pN N : ℝ) / dN N) z (E N)
        = ((((pN N : ℝ) / dN N : ℝ) : ℂ) - ((c : ℝ) : ℂ)) * E N := quad_sub_quad _ _ _ _
    have hcoe : ‖(((pN N : ℝ) / dN N : ℝ) : ℂ) - ((c : ℝ) : ℂ)‖
        = |(pN N : ℝ) / dN N - c| := by
      rw [← Complex.ofReal_sub, Complex.norm_real, Real.norm_eq_abs]
    have h2 : ‖MP.quad c z (E N)‖
        ≤ ‖MP.quad ((pN N : ℝ) / dN N) z (E N)‖
          + |(pN N : ℝ) / dN N - c| * ‖E N‖ := by
      have e : MP.quad c z (E N)
          = MP.quad ((pN N : ℝ) / dN N) z (E N)
            + ((((pN N : ℝ) / dN N : ℝ) : ℂ) - ((c : ℝ) : ℂ)) * E N := by
        rw [← hsub]; ring
      rw [e]
      refine (norm_add_le _ _).trans ?_
      rw [norm_mul, hcoe]
    have h3 : |(pN N : ℝ) / dN N - c| * ‖E N‖
        ≤ |(pN N : ℝ) / dN N - c| * (1 / z.im) :=
      mul_le_mul_of_nonneg_left hEnorm (abs_nonneg _)
    have h4 : resConst z / (dN N : ℝ) = resConst z / (dN N : ℕ) := rfl
    change ‖MP.quad c z (E N)‖ ≤ G N
    rw [hGdef]
    linarith
  have hRto : Tendsto R atTop (𝓝 0) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hGto
      (Filter.Eventually.of_forall fun N => norm_nonneg _) hRle
  have hRsmall : ∀ᶠ N in atTop, R N < 1 / 2 :=
    hRto.eventually (gt_mem_nhds (by norm_num : (0 : ℝ) < 1 / 2))
  -- the final bound, valid past the two eventual thresholds
  have hkey : ∀ᶠ N in atTop, ‖E N - MP.mC c z‖ ≤ R N / (‖z‖ * (z.im * rootLb c z ^ 2)) := by
    filter_upwards [hev1, hRsmall] with N hdN hRN
    set rr : ℂ := -(MP.quad c z (E N)) with hrr
    have hq : MP.quad c z (E N) = -rr := by rw [hrr, neg_neg]
    have hrrn : ‖rr‖ = R N := by rw [hrr, norm_neg]
    have hlb : rootLb c z ≤ ‖E N‖ :=
      norm_ge_rootLb hz0 hq (by rw [hrrn]; linarith)
    have him : z.im * ‖E N‖ ^ 2 ≤ (E N).im := im_integral_stieltjesC_ge hzim hdN
    have hsq : rootLb c z ^ 2 ≤ ‖E N‖ ^ 2 :=
      pow_le_pow_left₀ hδ.le hlb 2
    have h5 : z.im * rootLb c z ^ 2 ≤ z.im * ‖E N‖ ^ 2 :=
      mul_le_mul_of_nonneg_left hsq hzim.le
    have h6 : 0 < z.im * rootLb c z ^ 2 := by positivity
    have himpos : 0 < (E N).im := by linarith
    have hbnd := norm_sub_root_le hcpos hzim (MP.im_mC_pos hcpos.le hzim)
      (MP.quad_mC hcpos.le hzim) himpos hq
    rw [hrrn] at hbnd
    refine hbnd.trans ?_
    have hden1 : 0 < ‖z‖ * (E N).im := by positivity
    refine div_le_div_of_nonneg_left (norm_nonneg _) hDen ?_
    have hle : z.im * rootLb c z ^ 2 ≤ (E N).im := by linarith
    exact mul_le_mul_of_nonneg_left hle hzn.le
  have hdiv : Tendsto (fun N => R N / (‖z‖ * (z.im * rootLb c z ^ 2))) atTop (𝓝 0) := by
    simpa using hRto.div_const (‖z‖ * (z.im * rootLb c z ^ 2))
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hdiv
    (Filter.Eventually.of_forall fun N => norm_nonneg _) ?_
  simpa using hkey

/-! ### R1a and R1b -/

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  [∀ N, IsProbabilityMeasure (μ N)] {pN dN : ℕ → ℕ}

variable (hc : 0 < c) (hz : 0 < z.im) (hd : Tendsto dN atTop atTop) (hp : ∀ N, 0 < pN N)
  (hcN : Tendsto (fun N => (pN N : ℝ) / dN N) atTop (𝓝 c))
  (Y : ∀ N, Ω N → Matrix (Fin (pN N)) (Fin (dN N)) ℝ)
  (hY : ∀ N, HasLaw (Y N) (gaussianMatrix (pN N) (dN N)) (μ N))

omit [∀ N, IsProbabilityMeasure (μ N)] in
include hc hz hd hp hcN hY in
/-- **R1a.** `s_N(z) → MP.mC c z` in probability, in the norm form (choices 27, 28).

Proof: step 4 in the limit (`tendsto_norm_integral_sub_mC`) puts the mean within `ε/3` of the
root past some `N₀`; the two real deviations are each within `ε/3` outside a set of measure at
most `2 exp(-(ε/3)²/(2 L(z, d_N)²))` by `measure_abs_ge_le_of_lipschitz`; and `L(z, d_N) → 0`
because `d_N → ∞`. -/
theorem tendstoInProb_stieltjesC :
    TendstoInProb μ (fun N ω => ‖R4C.stieltjesC (W0 (Y N ω)) z - MP.mC c z‖) 0 := by
  intro ε hε
  have hz0 : z ≠ 0 := MP.ne_zero_of_im_pos hz
  set E : ℕ → ℂ :=
    fun N => ∫ Yv, R4C.stieltjesC (W0 Yv) z ∂gaussianMatrix (pN N) (dN N) with hEdef
  have hmean :=
    tendsto_norm_integral_sub_mC (c := c) (z := z) (pN := pN) (dN := dN) hc hz hd hp hcN
  have hev1 : ∀ᶠ N in atTop, 0 < dN N := hd.eventually_gt_atTop 0
  have hev2 : ∀ᶠ N in atTop, ‖E N - MP.mC c z‖ < ε / 3 :=
    hmean.eventually (gt_mem_nhds (by linarith : (0 : ℝ) < ε / 3))
  set K : ℝ := 2 * Real.sqrt (z.im + ‖z‖) / z.im ^ 2 with hKdef
  have hKpos : 0 < K := by
    have h1 : 0 < z.im + ‖z‖ := by linarith [norm_nonneg z]
    have h2 : 0 < Real.sqrt (z.im + ‖z‖) := Real.sqrt_pos.mpr h1
    rw [hKdef]; positivity
  set B : ℕ → ℝ≥0∞ :=
    fun N => ENNReal.ofReal (4 * Real.exp (-(ε / 3) ^ 2 / (2 * lipConst z (dN N) ^ 2)))
    with hBdef
  have hexpto : Tendsto (fun N => Real.exp (-(ε / 3) ^ 2 / (2 * lipConst z (dN N) ^ 2)))
      atTop (𝓝 0) := by
    refine Real.tendsto_exp_atBot.comp ?_
    have heq : (fun N => -(ε / 3) ^ 2 / (2 * lipConst z (dN N) ^ 2))
        =ᶠ[atTop] (fun N => (-(ε / 3) ^ 2 / (2 * K ^ 2)) * ((dN N : ℝ)) ^ 2) := by
      filter_upwards [hev1] with N hdN
      have hdR : (0 : ℝ) < dN N := by exact_mod_cast hdN
      have hl : lipConst z (dN N) = K / (dN N : ℝ) := by
        rw [hKdef]; unfold lipConst; ring
      rw [hl]
      field_simp
    refine Tendsto.congr' heq.symm ?_
    have hd2 : Tendsto (fun N => ((dN N : ℝ)) ^ 2) atTop atTop :=
      (tendsto_pow_atTop (α := ℝ) (n := 2) (by norm_num)).comp
        (tendsto_natCast_atTop_atTop.comp hd)
    have hneg : -(ε / 3) ^ 2 / (2 * K ^ 2) < 0 := by
      have h1 : (0 : ℝ) < (ε / 3) ^ 2 := by positivity
      have h2 : (0 : ℝ) < 2 * K ^ 2 := by positivity
      exact div_neg_of_neg_of_pos (by linarith) h2
    exact Filter.Tendsto.const_mul_atTop_of_neg hneg hd2
  have hBto : Tendsto B atTop (𝓝 0) := by
    have h4 : Tendsto (fun N => 4 * Real.exp (-(ε / 3) ^ 2 / (2 * lipConst z (dN N) ^ 2)))
        atTop (𝓝 0) := by simpa using hexpto.const_mul 4
    have := ENNReal.tendsto_ofReal h4
    simpa [hBdef] using this
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hBto
    (Filter.Eventually.of_forall fun N => zero_le) ?_
  filter_upwards [hev1, hev2] with N hdN hmeanN
  have hpd : 0 < pN N * dN N := Nat.mul_pos (hp N) hdN
  have hLpos : (0 : ℝ≥0) < lipConstNN z (dN N) :=
    NNReal.coe_pos.mp (lipConst_pos hz hdN)
  have hInt := integrable_stieltjesC (p := pN N) (d := dN N) hz hdN
  have hε3 : (0 : ℝ) < ε / 3 := by linarith
  -- the mean of the real and imaginary parts
  have hEre : (∫ Zv, (R4C.stieltjesC (W0 Zv) z).re ∂gaussianMatrix (pN N) (dN N)) = (E N).re := by
    have h := ContinuousLinearMap.integral_comp_comm Complex.reCLM hInt
    simpa [hEdef] using h
  have hEim : (∫ Zv, (R4C.stieltjesC (W0 Zv) z).im ∂gaussianMatrix (pN N) (dN N)) = (E N).im := by
    have h := ContinuousLinearMap.integral_comp_comm Complex.imCLM hInt
    simpa [hEdef] using h
  have hmre : MeasurableSet
      {Yv : Matrix (Fin (pN N)) (Fin (dN N)) ℝ |
        ε / 3 ≤ |(R4C.stieltjesC (W0 Yv) z).re - (E N).re|} :=
    measurableSet_le measurable_const
      (((measurable_stieltjesC_re hz hdN).sub measurable_const).abs)
  have hmim : MeasurableSet
      {Yv : Matrix (Fin (pN N)) (Fin (dN N)) ℝ |
        ε / 3 ≤ |(R4C.stieltjesC (W0 Yv) z).im - (E N).im|} :=
    measurableSet_le measurable_const
      (((measurable_stieltjesC_im hz hdN).sub measurable_const).abs)
  -- the two tail bounds on the canonical space
  have htre : gaussianMatrix (pN N) (dN N)
        {Yv | ε / 3 ≤ |(R4C.stieltjesC (W0 Yv) z).re - (E N).re|}
      ≤ ENNReal.ofReal (2 * Real.exp (-(ε / 3) ^ 2 / (2 * lipConst z (dN N) ^ 2))) := by
    have hkey := measure_abs_ge_le_of_lipschitz (p := pN N) (d := dN N) hpd hLpos
      (lipschitzWith_stieltjesC_re (p := pN N) (d := dN N) hz hdN) hε3
    have hkey2 : ((gaussianMatrix (pN N) (dN N))
        {Yv | ε / 3 ≤ |(R4C.stieltjesC (W0 Yv) z).re - (E N).re|}).toReal
          ≤ 2 * Real.exp (-(ε / 3) ^ 2 / (2 * lipConst z (dN N) ^ 2)) := by
      simpa [hEre] using hkey
    rw [← ENNReal.ofReal_toReal (measure_ne_top (gaussianMatrix (pN N) (dN N)) _)]
    exact ENNReal.ofReal_le_ofReal hkey2
  have htim : gaussianMatrix (pN N) (dN N)
        {Yv | ε / 3 ≤ |(R4C.stieltjesC (W0 Yv) z).im - (E N).im|}
      ≤ ENNReal.ofReal (2 * Real.exp (-(ε / 3) ^ 2 / (2 * lipConst z (dN N) ^ 2))) := by
    have hkey := measure_abs_ge_le_of_lipschitz (p := pN N) (d := dN N) hpd hLpos
      (lipschitzWith_stieltjesC_im (p := pN N) (d := dN N) hz hdN) hε3
    have hkey2 : ((gaussianMatrix (pN N) (dN N))
        {Yv | ε / 3 ≤ |(R4C.stieltjesC (W0 Yv) z).im - (E N).im|}).toReal
          ≤ 2 * Real.exp (-(ε / 3) ^ 2 / (2 * lipConst z (dN N) ^ 2)) := by
      simpa [hEim] using hkey
    rw [← ENNReal.ofReal_toReal (measure_ne_top (gaussianMatrix (pN N) (dN N)) _)]
    exact ENNReal.ofReal_le_ofReal hkey2
  -- the inclusion
  have hsub : {ω | ε ≤ |‖R4C.stieltjesC (W0 (Y N ω)) z - MP.mC c z‖ - 0|}
      ⊆ {ω | ε / 3 ≤ |(R4C.stieltjesC (W0 (Y N ω)) z).re - (E N).re|}
        ∪ {ω | ε / 3 ≤ |(R4C.stieltjesC (W0 (Y N ω)) z).im - (E N).im|} := by
    intro ω hω
    by_contra hcon
    simp only [Set.mem_union, Set.mem_ofPred_eq, not_or, not_le] at hcon
    obtain ⟨h1, h2⟩ := hcon
    have hmem : ε ≤ ‖R4C.stieltjesC (W0 (Y N ω)) z - MP.mC c z‖ := by
      have hω' : ε ≤ |‖R4C.stieltjesC (W0 (Y N ω)) z - MP.mC c z‖ - 0| := hω
      rwa [sub_zero, abs_of_nonneg (norm_nonneg _)] at hω'
    have hsplit : ‖R4C.stieltjesC (W0 (Y N ω)) z - E N‖
        ≤ |(R4C.stieltjesC (W0 (Y N ω)) z).re - (E N).re|
          + |(R4C.stieltjesC (W0 (Y N ω)) z).im - (E N).im| := by
      have h := Complex.norm_le_abs_re_add_abs_im (R4C.stieltjesC (W0 (Y N ω)) z - E N)
      simpa using h
    have htri : ‖R4C.stieltjesC (W0 (Y N ω)) z - MP.mC c z‖
        ≤ ‖R4C.stieltjesC (W0 (Y N ω)) z - E N‖ + ‖E N - MP.mC c z‖ := by
      have h := norm_add_le (R4C.stieltjesC (W0 (Y N ω)) z - E N) (E N - MP.mC c z)
      simpa using h
    linarith
  calc μ N {ω | ε ≤ |‖R4C.stieltjesC (W0 (Y N ω)) z - MP.mC c z‖ - 0|}
      ≤ μ N ({ω | ε / 3 ≤ |(R4C.stieltjesC (W0 (Y N ω)) z).re - (E N).re|}
          ∪ {ω | ε / 3 ≤ |(R4C.stieltjesC (W0 (Y N ω)) z).im - (E N).im|}) :=
        measure_mono hsub
    _ ≤ μ N {ω | ε / 3 ≤ |(R4C.stieltjesC (W0 (Y N ω)) z).re - (E N).re|}
          + μ N {ω | ε / 3 ≤ |(R4C.stieltjesC (W0 (Y N ω)) z).im - (E N).im|} :=
        measure_union_le _ _
    _ = gaussianMatrix (pN N) (dN N)
            {Yv | ε / 3 ≤ |(R4C.stieltjesC (W0 Yv) z).re - (E N).re|}
          + gaussianMatrix (pN N) (dN N)
            {Yv | ε / 3 ≤ |(R4C.stieltjesC (W0 Yv) z).im - (E N).im|} := by
        rw [(hY N).measure_eq hmre, (hY N).measure_eq hmim]
    _ ≤ ENNReal.ofReal (2 * Real.exp (-(ε / 3) ^ 2 / (2 * lipConst z (dN N) ^ 2)))
          + ENNReal.ofReal (2 * Real.exp (-(ε / 3) ^ 2 / (2 * lipConst z (dN N) ^ 2))) :=
        add_le_add htre htim
    _ = B N := by
        rw [hBdef, ← ENNReal.ofReal_add (by positivity) (by positivity)]
        ring_nf

omit [∀ N, IsProbabilityMeasure (μ N)] in
include hc hz hd hp hcN hY in
/-- **R1b.** `d⁻¹ tr G² → MP.mCDeriv c z`, same norm form.

Proof: Cauchy's estimate at radius `η/2` (`norm_stieltjes2C_sub_mCDeriv_le`) turns a uniform
bound on the circle into a bound on the derivative gap; the circle is compact, so a finite net
of radius `ρ` covers it; the gap is `32/η²`-Lipschitz in `ζ` on the disc of radius `3η/4`
(`norm_diff_sub_diff_le`), so the net values control the sup; and R1a at each of the finitely
many net points, with a union bound, sends the probability to `0`. -/
theorem tendstoInProb_stieltjes2C :
    TendstoInProb μ (fun N ω => ‖R4C.stieltjes2C (W0 (Y N ω)) z - MP.mCDeriv c z‖) 0 := by
  intro ε hε
  have hz0 : z ≠ 0 := MP.ne_zero_of_im_pos hz
  set r : ℝ := z.im / 2 with hrdef
  have hrpos : (0 : ℝ) < r := by rw [hrdef]; linarith
  set Lz : ℝ := 32 / z.im ^ 2 with hLzdef
  have hLzpos : (0 : ℝ) < Lz := by rw [hLzdef]; positivity
  set Csup : ℝ := ε * r / 2 with hCdef
  have hCpos : (0 : ℝ) < Csup := by rw [hCdef]; positivity
  set ρ : ℝ := Csup / (2 * Lz) with hrhodef
  have hrhopos : (0 : ℝ) < ρ := by rw [hrhodef]; positivity
  obtain ⟨b, hbsub, hbfin, hbcov⟩ :=
    (isCompact_sphere z r).elim_finite_subcover_image
      (b := Metric.sphere z r) (c := fun ζ : ℂ => Metric.ball ζ ρ)
      (fun ζ _ => Metric.isOpen_ball)
      (fun ζ hζ => Set.mem_biUnion hζ (Metric.mem_ball_self hrhopos))
  have hyim : ∀ y ∈ b, 0 < y.im := by
    intro y hy
    have hy' : y ∈ Metric.closedBall z r :=
      Metric.sphere_subset_closedBall (hbsub hy)
    have := im_ge_of_mem_closedBall (z := z) hy'
    rw [hrdef] at this
    linarith
  have hyball : ∀ y ∈ b, y ∈ Metric.closedBall z (3 * z.im / 4) := by
    intro y hy
    have h1 : dist y z = r := hbsub hy
    rw [Metric.mem_closedBall, h1, hrdef]
    linarith
  -- R1a at each net point
  have hpt : ∀ y ∈ hbfin.toFinset,
      Tendsto (fun N => μ N {ω | Csup / 2 ≤ ‖R4C.stieltjesC (W0 (Y N ω)) y - MP.mC c y‖})
        atTop (𝓝 0) := by
    intro y hy
    have hyb : y ∈ b := hbfin.mem_toFinset.mp hy
    have hR1 := tendstoInProb_stieltjesC (z := y) hc (hyim y hyb) hd hp hcN Y hY
      (Csup / 2) (by positivity)
    have hset : ∀ N, {ω : Ω N | Csup / 2 ≤ ‖R4C.stieltjesC (W0 (Y N ω)) y - MP.mC c y‖}
        = {ω | Csup / 2 ≤ |‖R4C.stieltjesC (W0 (Y N ω)) y - MP.mC c y‖ - 0|} := by
      intro N
      ext ω
      simp [abs_of_nonneg (norm_nonneg _)]
    simpa [hset] using hR1
  have hsum : Tendsto (fun N => ∑ y ∈ hbfin.toFinset,
      μ N {ω | Csup / 2 ≤ ‖R4C.stieltjesC (W0 (Y N ω)) y - MP.mC c y‖}) atTop (𝓝 0) := by
    have h := tendsto_finsetSum hbfin.toFinset hpt
    simpa using h
  have hev1 : ∀ᶠ N in atTop, 0 < dN N := hd.eventually_gt_atTop 0
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hsum
    (Filter.Eventually.of_forall fun N => zero_le) ?_
  filter_upwards [hev1] with N hdN
  have hsubset : {ω | ε ≤ |‖R4C.stieltjes2C (W0 (Y N ω)) z - MP.mCDeriv c z‖ - 0|}
      ⊆ ⋃ y ∈ hbfin.toFinset,
          {ω | Csup / 2 ≤ ‖R4C.stieltjesC (W0 (Y N ω)) y - MP.mC c y‖} := by
    intro ω hω
    by_contra hcon
    simp only [Set.mem_iUnion, Set.mem_ofPred_eq, not_exists, not_le] at hcon
    have hW := isHermitian_W0 (Y N ω)
    have hC : ∀ ζ ∈ Metric.sphere z r,
        ‖R4C.stieltjesC (W0 (Y N ω)) ζ - MP.mC c ζ‖ ≤ Csup := by
      intro ζ hζ
      obtain ⟨y, hyb, hyball'⟩ : ∃ y ∈ b, ζ ∈ Metric.ball y ρ := by
        have := hbcov hζ
        simpa using this
      have hlt : ‖R4C.stieltjesC (W0 (Y N ω)) y - MP.mC c y‖ < Csup / 2 :=
        hcon y (hbfin.mem_toFinset.mpr hyb)
      have hζball : ζ ∈ Metric.closedBall z (3 * z.im / 4) := by
        have h1 : dist ζ z = r := hζ
        rw [Metric.mem_closedBall, h1, hrdef]
        linarith
      have hlip := norm_diff_sub_diff_le (W := W0 (Y N ω)) hc hz hdN hW hζball
        (hyball y hyb)
      have hdistlt : ‖ζ - y‖ < ρ := by
        have : dist ζ y < ρ := hyball'
        rwa [Complex.dist_eq] at this
      have htri : ‖R4C.stieltjesC (W0 (Y N ω)) ζ - MP.mC c ζ‖
          ≤ ‖(R4C.stieltjesC (W0 (Y N ω)) ζ - MP.mC c ζ)
              - (R4C.stieltjesC (W0 (Y N ω)) y - MP.mC c y)‖
            + ‖R4C.stieltjesC (W0 (Y N ω)) y - MP.mC c y‖ := by
        have h := norm_add_le
          ((R4C.stieltjesC (W0 (Y N ω)) ζ - MP.mC c ζ)
            - (R4C.stieltjesC (W0 (Y N ω)) y - MP.mC c y))
          (R4C.stieltjesC (W0 (Y N ω)) y - MP.mC c y)
        simpa using h
      have hLzrho : Lz * ρ = Csup / 2 := by
        rw [hrhodef]
        field_simp
      have hstep : Lz * ‖ζ - y‖ ≤ Lz * ρ :=
        mul_le_mul_of_nonneg_left hdistlt.le hLzpos.le
      have hlip' : ‖(R4C.stieltjesC (W0 (Y N ω)) ζ - MP.mC c ζ)
          - (R4C.stieltjesC (W0 (Y N ω)) y - MP.mC c y)‖ ≤ Lz * ‖ζ - y‖ := by
        rw [hLzdef]; exact hlip
      linarith
    have hcauchy := norm_stieltjes2C_sub_mCDeriv_le (W := W0 (Y N ω)) hc hz hdN hW hC
    have hmem : ε ≤ ‖R4C.stieltjes2C (W0 (Y N ω)) z - MP.mCDeriv c z‖ := by
      have hω' : ε ≤ |‖R4C.stieltjes2C (W0 (Y N ω)) z - MP.mCDeriv c z‖ - 0| := hω
      rwa [sub_zero, abs_of_nonneg (norm_nonneg _)] at hω'
    have hval : Csup / (z.im / 2) = ε / 2 := by
      rw [hCdef, hrdef]
      field_simp
    rw [← hrdef, hval] at hcauchy
    linarith
  exact (measure_mono hsubset).trans (measure_biUnion_finset_le _ _)

end R1

end StackedSVD


