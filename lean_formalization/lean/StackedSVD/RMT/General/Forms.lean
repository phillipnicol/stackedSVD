/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.ProbC
import StackedSVD.RMT.General.Stability
import StackedSVD.RMT.General.FormsBridge

/-!
# The derivative step for the bilinear resolvent forms, general noise law

Stage 1, unit F1a, part 2. `RMT/General/Trace.lean` sends `d⁻¹ tr G(z)` to `m_c(z)` and then
`d⁻¹ tr G(z)²` to `m_c'(z)` by a Cauchy estimate. The same step is needed for the bilinear
forms `Φ_v`, `Φ_u` and `Ψ_{v,g}`. This file proves it once, in a model free form, and applies
it to the three forms.

## Content

1. `hasDerivAt_cformC`: `ζ ↦ x ⬝ᵥ G(ζ) y` has derivative `x ⬝ᵥ G(z)² y` at every `z` above
   the real axis. The eigen-sum of `R4C.cformC_eq_sum` is differentiated term by term, and
   `R4C.cform2C_eq_sum` names the sum of the derivatives. The power convention matches: the
   derivative of `(λ - ζ)⁻¹` is `((λ - z)⁻¹)²`, with no sign and no extra factor, which is
   the summand of `cform2C_eq_sum`.
2. `tendstoInProbC_deriv_of_forms`: the generic transfer. If a family of functions `f N ω`
   holomorphic above the real axis converges in probability to `F` at every point above the
   axis, and the gap of the derivatives is bounded by one constant `L` on a disc around `z`
   with probability tending to one, then the derivatives converge in probability at `z`. The
   route is the one of `tendstoInProb_stieltjes2C_general` (`RMT/General/Trace.lean:1133`):
   the mean value inequality on the convex disc, Cauchy's estimate on the sphere of radius
   `Im z / 2`, a finite net of the (compact) sphere, and a union bound.
3. The three instances, in `namespace SpikedModel`: `qform2C_gram_v_tendsto`,
   `qform2C_gramC_u_tendsto` and `cform2C_gram_v_gvec_tendsto`. Each takes the matching
   first order limit as a hypothesis and returns the second order one at the same point.

Choice 8 of `notes/archive/prop_single_table_general.md` keeps every file under `RMT/General/` free
of `RMT/R1.lean`, `RMT/R2.lean`, `RMT/SteinStep.lean` and `Vendor/COLT83/`; this file follows
that boundary and uses the `GenRMT` copies of `RMT/General/Stability.lean`.
-/

open MeasureTheory Filter Topology Set
open scoped Matrix

namespace StackedSVD
namespace GenRMT

/-! ### F1: the derivative of the bilinear resolvent form -/

variable {d : ℕ} {W : Matrix (Fin d) (Fin d) ℝ} {z : ℂ}

/-- **The derivative of `Ψ_{x,y}` in `z` is `Ψ²_{x,y}`.** The bilinear twin of
`hasDerivAt_stieltjesC` (`RMT/General/Stability.lean:98`, itself a copy of `R1.lean:227`).
The convention is the literal one: `R4C.cform2C_eq_sum` sums `((λ_a - z)⁻¹)²` against the
coefficient `(Uᵀx)_a (Uᵀy)_a`, and that is exactly the term by term derivative of the sum of
`R4C.cformC_eq_sum`, with no sign and no extra factor. -/
theorem hasDerivAt_cformC (hW : W.IsHermitian) (hz : 0 < z.im) (x y : Fin d → ℝ) :
    HasDerivAt (fun ζ => R4C.cformC W ζ x y) (R4C.cform2C W z x y) z := by
  have hopen : IsOpen {ζ : ℂ | 0 < ζ.im} := isOpen_lt continuous_const Complex.continuous_im
  have hmem : z ∈ {ζ : ℂ | 0 < ζ.im} := hz
  have hkey : HasDerivAt (fun ζ : ℂ => ∑ i, ((hW.eigenvalues i : ℂ) - ζ)⁻¹ *
        ((((R4.eigU hW)ᵀ *ᵥ x) i * ((R4.eigU hW)ᵀ *ᵥ y) i : ℝ) : ℂ))
      (∑ i, (((hW.eigenvalues i : ℂ) - z)⁻¹) ^ 2 *
        ((((R4.eigU hW)ᵀ *ᵥ x) i * ((R4.eigU hW)ᵀ *ᵥ y) i : ℝ) : ℂ)) z := by
    refine HasDerivAt.fun_sum fun i _ => ?_
    have h1 : HasDerivAt (fun ζ : ℂ => (hW.eigenvalues i : ℂ) - ζ) (-1) z :=
      HasDerivAt.const_sub ((hW.eigenvalues i : ℂ)) (hasDerivAt_id z)
    have h2 := h1.inv (R4C.eigenvalue_sub_ne_zero hW hz.ne' i)
    have heq : -(-1 : ℂ) / ((hW.eigenvalues i : ℂ) - z) ^ 2
        = (((hW.eigenvalues i : ℂ) - z)⁻¹) ^ 2 := by
      rw [inv_pow, neg_neg, one_div]
    rw [← heq]
    exact h2.mul_const _
  rw [R4C.cform2C_eq_sum hW hz.ne']
  refine hkey.congr_of_eventuallyEq ?_
  filter_upwards [hopen.mem_nhds hmem] with ζ hζ
  show R4C.cformC W ζ x y = _
  rw [R4C.cformC_eq_sum hW (ne_of_gt hζ)]

/-- The derivative of `Φ_y` in `z` is `Φ²_y`, the diagonal case of `hasDerivAt_cformC`. -/
theorem hasDerivAt_qformC (hW : W.IsHermitian) (hz : 0 < z.im) (y : Fin d → ℝ) :
    HasDerivAt (fun ζ => R4C.qformC W ζ y) (R4C.qform2C W z y) z :=
  hasDerivAt_cformC hW hz y y

/-! ### F2: the generic transfer from a limit of forms to a limit of derivatives -/

section Transfer

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}

/-- **The Cauchy transfer, model free.** Let `f N ω` be holomorphic above the real axis with
derivative `f' N ω`, and `F` holomorphic there with derivative `F'`. Assume:

* `f N ω ζ → F ζ` in probability at every `ζ` above the axis;
* on a good event `S N` of probability tending to one, the gap of the derivatives is at most
  one constant `L` on the closed disc of radius `3 Im z / 4` around `z`.

Then the gap of the derivatives at `z` tends to `0` in probability.

The route is the one of `tendstoInProb_stieltjes2C_general` (`RMT/General/Trace.lean:1133`),
with the model specific pieces replaced by the two hypotheses: the difference is Lipschitz on
the convex disc by the mean value inequality, so a finite `ρ`-net of the compact sphere of
radius `Im z / 2` turns finitely many pointwise limits into a uniform bound on the sphere, and
Cauchy's estimate at that radius turns the uniform bound into a bound on the derivative gap at
the center. `L` is only used as an upper bound, so a value that is not positive is harmless:
the proof runs with `max L 1`. -/
theorem tendstoInProbC_deriv_of_forms (f f' : ∀ N, Ω N → ℂ → ℂ) (F F' : ℂ → ℂ)
    {z : ℂ} (hz : 0 < z.im) (L : ℝ) (S : ∀ N, Set (Ω N))
    (hS : Tendsto (fun N => μ N (S N)ᶜ) atTop (𝓝 0))
    (hderiv : ∀ N ω ζ, 0 < ζ.im → HasDerivAt (f N ω) (f' N ω ζ) ζ)
    (hF : ∀ ζ, 0 < ζ.im → HasDerivAt F (F' ζ) ζ)
    (hbound : ∀ N ω, ω ∈ S N → ∀ ζ ∈ Metric.closedBall z (3 * z.im / 4),
      ‖f' N ω ζ - F' ζ‖ ≤ L)
    (hlim : ∀ ζ, 0 < ζ.im → TendstoInProb μ (fun N ω => ‖f N ω ζ - F ζ‖) 0) :
    TendstoInProb μ (fun N ω => ‖f' N ω z - F' z‖) 0 := by
  classical
  intro ε hε
  set r : ℝ := z.im / 2 with hrdef
  have hrpos : (0 : ℝ) < r := by rw [hrdef]; linarith
  set Lz : ℝ := max L 1 with hLzdef
  have hLzpos : (0 : ℝ) < Lz := lt_of_lt_of_le one_pos (le_max_right _ _)
  set Csup : ℝ := ε * r / 2 with hCdef
  have hCpos : (0 : ℝ) < Csup := by rw [hCdef]; positivity
  set ρ : ℝ := Csup / (2 * Lz) with hrhodef
  have hrhopos : (0 : ℝ) < ρ := by rw [hrhodef]; positivity
  obtain ⟨b, hbsub, hbfin, hbcov⟩ :=
    (isCompact_sphere z r).elim_finite_subcover_image
      (b := Metric.sphere z r) (c := fun ζ : ℂ => Metric.ball ζ ρ)
      (fun ζ _ => Metric.isOpen_ball)
      (fun ζ hζ => Set.mem_biUnion hζ (Metric.mem_ball_self hrhopos))
  have hyball : ∀ y ∈ b, y ∈ Metric.closedBall z (3 * z.im / 4) := by
    intro y hy
    have h1 : dist y z = r := hbsub hy
    rw [Metric.mem_closedBall, h1, hrdef]
    linarith
  have hyim : ∀ y ∈ b, 0 < y.im := by
    intro y hy
    have := im_ge_of_mem_closedBall (z := z) (hyball y hy)
    linarith
  have hpt : ∀ y ∈ hbfin.toFinset,
      Tendsto (fun N => μ N {ω | Csup / 2 ≤ ‖f N ω y - F y‖}) atTop (𝓝 0) := fun y hy =>
    tendstoInProbC_iff.mp (hlim y (hyim y (hbfin.mem_toFinset.mp hy))) (Csup / 2) (by positivity)
  have hsum : Tendsto (fun N => μ N (S N)ᶜ +
      ∑ y ∈ hbfin.toFinset, μ N {ω | Csup / 2 ≤ ‖f N ω y - F y‖}) atTop (𝓝 0) := by
    have hfin : Tendsto (fun N => ∑ y ∈ hbfin.toFinset,
        μ N {ω | Csup / 2 ≤ ‖f N ω y - F y‖}) atTop (𝓝 0) := by
      simpa using tendsto_finsetSum hbfin.toFinset hpt
    simpa using hS.add hfin
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hsum
    (fun _ => zero_le) fun N => ?_
  have hsubset : {ω | ε ≤ |‖f' N ω z - F' z‖ - 0|}
      ⊆ (S N)ᶜ ∪ ⋃ y ∈ hbfin.toFinset, {ω | Csup / 2 ≤ ‖f N ω y - F y‖} := by
    intro ω hω
    by_cases hωS : ω ∈ S N
    · refine Set.mem_union_right _ ?_
      by_contra hcon
      simp only [Set.mem_iUnion, Set.mem_ofPred_eq, not_exists, not_le] at hcon
      -- the difference and its derivative
      set D : ℂ → ℂ := fun w => f N ω w - F w with hDdef
      have hDder : ∀ w : ℂ, 0 < w.im → HasDerivAt D (f' N ω w - F' w) w := fun w hw =>
        (hderiv N ω w hw).sub (hF w hw)
      have himball : ∀ w ∈ Metric.closedBall z (3 * z.im / 4), z.im / 4 ≤ w.im := by
        intro w hw
        have := im_ge_of_mem_closedBall (z := z) hw
        linarith
      have hbnd : ∀ w ∈ Metric.closedBall z (3 * z.im / 4), ‖f' N ω w - F' w‖ ≤ Lz :=
        fun w hw => (hbound N ω hωS w hw).trans (le_max_left _ _)
      have hderW : ∀ w ∈ Metric.closedBall z (3 * z.im / 4),
          HasDerivWithinAt D (f' N ω w - F' w) (Metric.closedBall z (3 * z.im / 4)) w := by
        intro w hw
        have h1 : z.im / 4 ≤ w.im := himball w hw
        exact (hDder w (by linarith)).hasDerivWithinAt
      -- a uniform bound on the sphere of radius `r`
      have hC : ∀ ζ ∈ Metric.sphere z r, ‖D ζ‖ ≤ Csup := by
        intro ζ hζ
        obtain ⟨y, hyb, hyball'⟩ : ∃ y ∈ b, ζ ∈ Metric.ball y ρ := by
          have := hbcov hζ
          simpa using this
        have hlt : ‖f N ω y - F y‖ < Csup / 2 := by
          have := hcon y (hbfin.mem_toFinset.mpr hyb)
          simpa using this
        have hζball : ζ ∈ Metric.closedBall z (3 * z.im / 4) := by
          have h1 : dist ζ z = r := hζ
          rw [Metric.mem_closedBall, h1, hrdef]
          linarith
        have hcvx : Convex ℝ (Metric.closedBall z (3 * z.im / 4)) := convex_closedBall _ _
        have hlip := hcvx.norm_image_sub_le_of_norm_hasDerivWithin_le hderW hbnd
          (hyball y hyb) hζball
        have hdistlt : ‖ζ - y‖ < ρ := by
          have hd' : dist ζ y < ρ := hyball'
          rwa [Complex.dist_eq] at hd'
        have hLzrho : Lz * ρ = Csup / 2 := by
          rw [hrhodef]
          field_simp
        have hstep : Lz * ‖ζ - y‖ ≤ Lz * ρ := mul_le_mul_of_nonneg_left hdistlt.le hLzpos.le
        have htri : ‖D ζ‖ ≤ ‖D ζ - D y‖ + ‖D y‖ := by
          simpa using norm_add_le (D ζ - D y) (D y)
        have hDy : ‖D y‖ = ‖f N ω y - F y‖ := rfl
        rw [hDy] at htri
        linarith
      -- Cauchy's estimate at radius `r`
      have hdiffOn : DifferentiableOn ℂ D (Metric.closedBall z r) := by
        intro ζ hζ
        have h1 : ζ ∈ Metric.closedBall z (3 * z.im / 4) := by
          refine Metric.closedBall_subset_closedBall ?_ hζ
          rw [hrdef]; linarith
        exact ((hDder ζ (by linarith [himball ζ h1])).differentiableAt).differentiableWithinAt
      have hdc : DiffContOnCl ℂ D (Metric.ball z r) := by
        constructor
        · exact hdiffOn.mono Metric.ball_subset_closedBall
        · rw [closure_ball z (ne_of_gt hrpos)]
          exact hdiffOn.continuousOn
      have hderivz : deriv D z = f' N ω z - F' z := (hDder z hz).deriv
      have hcauchy : ‖f' N ω z - F' z‖ ≤ Csup / r := by
        rw [← hderivz]
        exact Complex.norm_deriv_le_of_forall_mem_sphere_norm_le hrpos hdc hC
      have hmem : ε ≤ ‖f' N ω z - F' z‖ := by
        have hω' : ε ≤ |‖f' N ω z - F' z‖ - 0| := hω
        rwa [sub_zero, abs_of_nonneg (norm_nonneg _)] at hω'
      have hval : Csup / r = ε / 2 := by
        rw [hCdef]
        field_simp
      rw [hval] at hcauchy
      linarith
    · exact Set.mem_union_left _ hωS
  calc μ N {ω | ε ≤ |‖f' N ω z - F' z‖ - 0|}
      ≤ μ N ((S N)ᶜ ∪ ⋃ y ∈ hbfin.toFinset, {ω | Csup / 2 ≤ ‖f N ω y - F y‖}) :=
        measure_mono hsubset
    _ ≤ μ N (S N)ᶜ + μ N (⋃ y ∈ hbfin.toFinset, {ω | Csup / 2 ≤ ‖f N ω y - F y‖}) :=
        measure_union_le _ _
    _ ≤ μ N (S N)ᶜ + ∑ y ∈ hbfin.toFinset, μ N {ω | Csup / 2 ≤ ‖f N ω y - F y‖} :=
        add_le_add le_rfl (measure_biUnion_finset_le _ _)

end Transfer

/-! ### A uniform bound on the disc

The two resolvent bounds and the two limit bounds are all of the form `1 / (Im ζ)²`. On the
disc that `tendstoInProbC_deriv_of_forms` uses, `Im ζ ≥ Im z / 4`, so each is at most
`16 / (Im z)²`. -/

/-- On the closed disc of radius `3 Im z / 4` around `z` the imaginary part stays positive and
`1 / (Im ζ)² ≤ 16 / (Im z)²`. -/
theorem one_div_sq_im_le_of_mem_closedBall (hz : 0 < z.im) {ζ : ℂ}
    (hζ : ζ ∈ Metric.closedBall z (3 * z.im / 4)) :
    0 < ζ.im ∧ 1 / ζ.im ^ 2 ≤ 16 / z.im ^ 2 := by
  have h1 : z.im / 4 ≤ ζ.im := by
    have := im_ge_of_mem_closedBall (z := z) hζ
    linarith
  have hpos : 0 < ζ.im := by linarith
  refine ⟨hpos, ?_⟩
  rw [div_le_div_iff₀ (by positivity) (by positivity)]
  have hkey : (0 : ℝ) ≤ (4 * ζ.im - z.im) * (4 * ζ.im + z.im) :=
    mul_nonneg (by linarith) (by linarith)
  nlinarith [hkey]

end GenRMT

/-! ### F3: the three instances -/

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- **L1'.** From `Φ_v(z) → m_c(z)` at every point above the axis to `Φ²_v(z) → m_c'(z)` at
one such point. The good event is everything: `‖Φ²_v(ζ)‖ ≤ (v ⬝ v) / (Im ζ)² = 1 / (Im ζ)²`
holds for every realization (`R4C.norm_qform2C_le`, `SpikedModel.dotProduct_v_self`), and
`‖m_c'(ζ)‖ ≤ 1 / (Im ζ)²` is `MP.norm_mCDeriv_le'`. -/
theorem qform2C_gram_v_tendsto {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (_hreg : m.Regime c) {ν : Measure ℝ} (_hν : NoiseLaw ν)
    (_hG : m.GeneralNoise ν)
    (hL : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω =>
      ‖R4C.qformC (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N)) - MP.mC c z‖) 0)
    (z : ℂ) (hz : 0 < z.im) :
    TendstoInProb μ (fun N ω =>
      ‖R4C.qform2C (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N)) - MP.mCDeriv c z‖) 0 := by
  refine GenRMT.tendstoInProbC_deriv_of_forms
    (fun N ω ζ => R4C.qformC (GenRMT.gram (m.Z N ω)) ζ (WithLp.ofLp (m.v N)))
    (fun N ω ζ => R4C.qform2C (GenRMT.gram (m.Z N ω)) ζ (WithLp.ofLp (m.v N)))
    (MP.mC c) (MP.mCDeriv c) hz (32 / z.im ^ 2) (fun _ => Set.univ) ?_ ?_ ?_ ?_ hL
  · simp
  · intro N ω ζ hζ
    exact GenRMT.hasDerivAt_qformC (GenRMT.gram_isHermitian _) hζ _
  · intro ζ hζ
    exact MP.hasDerivAt_mC hc.le hζ
  · intro N ω _ ζ hζ
    obtain ⟨hζpos, hζle⟩ := GenRMT.one_div_sq_im_le_of_mem_closedBall hz hζ
    have h1 : ‖R4C.qform2C (GenRMT.gram (m.Z N ω)) ζ (WithLp.ofLp (m.v N))‖
        ≤ 16 / z.im ^ 2 := by
      have hb := R4C.norm_qform2C_le (GenRMT.gram_isHermitian (m.Z N ω)) hζpos
        (WithLp.ofLp (m.v N))
      rw [m.dotProduct_v_self N] at hb
      exact hb.trans hζle
    have h2 : ‖MP.mCDeriv c ζ‖ ≤ 16 / z.im ^ 2 := (MP.norm_mCDeriv_le' hc hζpos).trans hζle
    calc ‖R4C.qform2C (GenRMT.gram (m.Z N ω)) ζ (WithLp.ofLp (m.v N)) - MP.mCDeriv c ζ‖
        ≤ ‖R4C.qform2C (GenRMT.gram (m.Z N ω)) ζ (WithLp.ofLp (m.v N))‖
            + ‖MP.mCDeriv c ζ‖ := norm_sub_le _ _
      _ ≤ 32 / z.im ^ 2 := by
          have hsplit : (32 : ℝ) / z.im ^ 2 = 16 / z.im ^ 2 + 16 / z.im ^ 2 := by ring
          rw [hsplit]
          exact add_le_add h1 h2

/-- **L2'.** The companion twin of L1': from `Φ_u(z) → m̃_c(z)` to `Φ²_u(z) → m̃_c'(z)`. The
bound on `m̃_c'` comes from `MP.mTildeCDeriv_eq_mCDeriv_inv` (`FormsBridge.lean:210`), which
reads `m̃_c'(ζ)` as `c⁻² m_{c⁻¹}'(ζ / c)`; with `(ζ / c).im = ζ.im / c` this is again at most
`1 / (Im ζ)²`. -/
theorem qform2C_gramC_u_tendsto {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (_hreg : m.Regime c) {ν : Measure ℝ} (_hν : NoiseLaw ν)
    (_hG : m.GeneralNoise ν)
    (hL : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω =>
      ‖R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)) - MP.mTildeC c z‖) 0)
    (z : ℂ) (hz : 0 < z.im) :
    TendstoInProb μ (fun N ω =>
      ‖R4C.qform2C (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))
        - MP.mTildeCDeriv c z‖) 0 := by
  refine GenRMT.tendstoInProbC_deriv_of_forms
    (fun N ω ζ => R4C.qformC (GenRMT.gramC (m.Z N ω)) ζ (WithLp.ofLp (m.u N)))
    (fun N ω ζ => R4C.qform2C (GenRMT.gramC (m.Z N ω)) ζ (WithLp.ofLp (m.u N)))
    (MP.mTildeC c) (MP.mTildeCDeriv c) hz (32 / z.im ^ 2) (fun _ => Set.univ) ?_ ?_ ?_ ?_ hL
  · simp
  · intro N ω ζ hζ
    exact GenRMT.hasDerivAt_qformC (GenRMT.gramC_isHermitian _) hζ _
  · intro ζ hζ
    exact MP.hasDerivAt_mTildeC hc hζ
  · intro N ω _ ζ hζ
    obtain ⟨hζpos, hζle⟩ := GenRMT.one_div_sq_im_le_of_mem_closedBall hz hζ
    have h1 : ‖R4C.qform2C (GenRMT.gramC (m.Z N ω)) ζ (WithLp.ofLp (m.u N))‖
        ≤ 16 / z.im ^ 2 := by
      have hb := R4C.norm_qform2C_le (GenRMT.gramC_isHermitian (m.Z N ω)) hζpos
        (WithLp.ofLp (m.u N))
      rw [m.dotProduct_u_self N] at hb
      exact hb.trans hζle
    have h2 : ‖MP.mTildeCDeriv c ζ‖ ≤ 16 / z.im ^ 2 := by
      have hzc : 0 < (ζ / (c : ℂ)).im := by
        rw [Complex.div_ofReal_im]
        exact div_pos hζpos hc
      have hb := MP.norm_mCDeriv_le' (c := c⁻¹) (inv_pos.mpr hc) hzc
      rw [Complex.div_ofReal_im] at hb
      have hnc : ‖(((c : ℝ) : ℂ)⁻¹) ^ 2‖ = (c⁻¹) ^ 2 := by
        rw [norm_pow, norm_inv, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hc]
      have hstep : ‖MP.mTildeCDeriv c ζ‖ ≤ 1 / ζ.im ^ 2 := by
        rw [MP.mTildeCDeriv_eq_mCDeriv_inv hc hζpos, norm_mul, hnc]
        calc (c⁻¹) ^ 2 * ‖MP.mCDeriv c⁻¹ (ζ / (c : ℂ))‖
            ≤ (c⁻¹) ^ 2 * (1 / (ζ.im / c) ^ 2) :=
              mul_le_mul_of_nonneg_left hb (by positivity)
          _ = 1 / ζ.im ^ 2 := by
              field_simp
      exact hstep.trans hζle
    calc ‖R4C.qform2C (GenRMT.gramC (m.Z N ω)) ζ (WithLp.ofLp (m.u N)) - MP.mTildeCDeriv c ζ‖
        ≤ ‖R4C.qform2C (GenRMT.gramC (m.Z N ω)) ζ (WithLp.ofLp (m.u N))‖
            + ‖MP.mTildeCDeriv c ζ‖ := norm_sub_le _ _
      _ ≤ 32 / z.im ^ 2 := by
          have hsplit : (32 : ℝ) / z.im ^ 2 = 16 / z.im ^ 2 + 16 / z.im ^ 2 := by ring
          rw [hsplit]
          exact add_le_add h1 h2

/-- **L3'.** The cross form: from `Ψ_{v,g}(z) → 0` to `Ψ²_{v,g}(z) → 0`. Here the limit and
its derivative are both `0`, and the good event is `g ⬝ g ≤ 4`, whose complement has
vanishing measure because `g ⬝ g → 1` in probability
(`tendstoInProb_dotProduct_gvec_general`, `RMT/General/Defs.lean:239`). On it
`‖Ψ²_{v,g}(ζ)‖ ≤ √(v ⬝ v) √(g ⬝ g) / (Im ζ)² ≤ 2 / (Im ζ)²`. -/
theorem cform2C_gram_v_gvec_tendsto {c : ℝ} (_hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hG : m.GeneralNoise ν)
    (hL : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω =>
      ‖R4C.cformC (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N)) (m.gvec N ω)‖) 0)
    (z : ℂ) (hz : 0 < z.im) :
    TendstoInProb μ (fun N ω =>
      ‖R4C.cform2C (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N)) (m.gvec N ω)‖) 0 := by
  have hd : Tendsto d atTop atTop := hreg.2.1
  have hgdot : TendstoInProb μ (fun N ω => m.gvec N ω ⬝ᵥ m.gvec N ω) 1 :=
    m.tendstoInProb_dotProduct_gvec_general hν hG hd
  have hgbad :
      Tendsto (fun N => μ N {ω : Ω N | m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4}ᶜ) atTop (𝓝 0) := by
    refine tendsto_measure_zero_of_subset
      (t := fun N => {ω : Ω N | (3 : ℝ) ≤ |m.gvec N ω ⬝ᵥ m.gvec N ω - 1|})
      (fun N ω hω => ?_) (hgdot 3 (by norm_num))
    have hnot : ¬ (m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4) := hω
    have h4 : (4 : ℝ) < m.gvec N ω ⬝ᵥ m.gvec N ω := not_le.mp hnot
    change (3 : ℝ) ≤ |m.gvec N ω ⬝ᵥ m.gvec N ω - 1|
    rw [abs_of_nonneg (by linarith)]
    linarith
  have hsq4 : Real.sqrt 4 = 2 := by
    rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num : (0 : ℝ) ≤ 2)]
  have key : TendstoInProb μ (fun N ω =>
      ‖R4C.cform2C (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N)) (m.gvec N ω)
        - (0 : ℂ)‖) 0 := by
    refine GenRMT.tendstoInProbC_deriv_of_forms
      (fun N ω ζ => R4C.cformC (GenRMT.gram (m.Z N ω)) ζ (WithLp.ofLp (m.v N)) (m.gvec N ω))
      (fun N ω ζ => R4C.cform2C (GenRMT.gram (m.Z N ω)) ζ (WithLp.ofLp (m.v N)) (m.gvec N ω))
      (fun _ => 0) (fun _ => 0) hz (32 / z.im ^ 2)
      (fun N => {ω : Ω N | m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4}) hgbad ?_ ?_ ?_ ?_
    · intro N ω ζ hζ
      exact GenRMT.hasDerivAt_cformC (GenRMT.gram_isHermitian _) hζ _ _
    · intro ζ _
      exact hasDerivAt_const ζ 0
    · intro N ω hωS ζ hζ
      obtain ⟨hζpos, hζle⟩ := GenRMT.one_div_sq_im_le_of_mem_closedBall hz hζ
      have hgood : m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4 := hωS
      have hb := R4C.norm_cform2C_le (GenRMT.gram_isHermitian (m.Z N ω)) hζpos
        (WithLp.ofLp (m.v N)) (m.gvec N ω)
      rw [m.dotProduct_v_self N, Real.sqrt_one, one_mul] at hb
      have hs : Real.sqrt (m.gvec N ω ⬝ᵥ m.gvec N ω) ≤ 2 := by
        rw [← hsq4]
        exact Real.sqrt_le_sqrt hgood
      have hstep : Real.sqrt (m.gvec N ω ⬝ᵥ m.gvec N ω) / ζ.im ^ 2 ≤ 32 / z.im ^ 2 := by
        have h2 : Real.sqrt (m.gvec N ω ⬝ᵥ m.gvec N ω) / ζ.im ^ 2
            ≤ 2 * (1 / ζ.im ^ 2) := by
          rw [div_eq_mul_one_div]
          exact mul_le_mul_of_nonneg_right hs (by positivity)
        have h3 : 2 * (1 / ζ.im ^ 2) ≤ 2 * (16 / z.im ^ 2) :=
          mul_le_mul_of_nonneg_left hζle (by norm_num)
        calc Real.sqrt (m.gvec N ω ⬝ᵥ m.gvec N ω) / ζ.im ^ 2
            ≤ 2 * (1 / ζ.im ^ 2) := h2
          _ ≤ 2 * (16 / z.im ^ 2) := h3
          _ = 32 / z.im ^ 2 := by ring
      rw [sub_zero]
      exact hb.trans hstep
    · intro ζ hζ
      simpa using hL ζ hζ
  simpa using key

end SpikedModel

end StackedSVD

