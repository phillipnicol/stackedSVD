/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Iso
import StackedSVD.RMT.General.Trace
import StackedSVD.RMT.General.FormsBridge
import StackedSVD.RMT.General.MPtilde
import StackedSVD.RMT.General.R3minus

/-!
# Unit F2: the two quadratic-form limits at a general noise law

`notes/archive/prop_single_table_general.md` section 5, unit F2 (Stage 1). The two scalar limits L1
and L2 of `RMT/General/FormsGeneral.lean`'s stub table: the quadratic form of the resolvent at
the deterministic directions `v` (the Wishart block) and `u` (the companion block) converges to
the Marchenko-Pastur transform and to its companion transform, in probability.

* L1, `qformC_gram_v_tendsto`: `vᵀ (gram Z - z)⁻¹ v → MP.mC c z`.
* L2, `qformC_gramC_u_tendsto`: `uᵀ (gramC Z - z)⁻¹ u → MP.mTildeC c z`.

## Route

The private lemma `qform_tendsto_aux` is the shared core of L1 and L2: at a general noise
matrix `Y` of shape `pN × dN` with i.i.d. entries of law `ν`, along a regime `pN / dN → c'`,
the quadratic form of the resolvent at a sequence of unit vectors converges to `MP.mC c' z` in
probability. It combines two facts already at a general law by the triangle inequality.

1. **Statement 1** of item G4 (`GenRMT.measure_qformC_sub_ge_le`, `RMT/General/Iso.lean`): at
   fixed `pN N`, `dN N`, the quadratic form `qformC (gram Y) z u` is within `ε` of
   `(u ⬝ u) stieltjesC (gram Y) z` outside a set of measure at most `isoRate ν z (pN N) (dN N) /
   ε²`; `u ⬝ u = 1` drops the coefficient. `GenRMT.tendsto_isoRate` sends the rate to `0` along
   the regime, so `Cheb.tendstoInProb_of_meas_le` turns the bound into a `TendstoInProb`
   statement after the law `hY` moves it from `noiseMatrix` to `μ N`.
2. **The trace law** of item G3 (`GenRMT.tendstoInProb_stieltjesC_general`,
   `RMT/General/Trace.lean`): `stieltjesC (gram Y) z → MP.mC c' z` in probability.
3. The triangle inequality `‖qformC - mC‖ ≤ ‖qformC - stieltjesC‖ + ‖stieltjesC - mC‖` and
   `TendstoInProb.add`/`TendstoInProb.of_le` (`Prob/TendstoInProb.lean`) close the lemma.

L1 is `qform_tendsto_aux` at `Y := m.Z`, `u := m.v`, `c' := c`.

L2 reduces to L1's core by duality: the companion form `qformC (gramC Z) z u` is a rescaled
form of `gram (Zᵀ)` at the point `z d / n` (`GenRMT.cformC_gramC_eq`, C2 of
`RMT/General/FormsBridge.lean`), `Zᵀ` again has i.i.d. entries of law `ν`
(`GenRMT.measurePreserving_transpose_noiseMatrix`), and `d / n → c⁻¹`, so `qform_tendsto_aux`
at the transposed model and the FIXED point `z' := z / c` gives `qformC (gram Zᵀ) z' u →
MP.mC c⁻¹ z'`; the Lipschitz bound of `cformC` in `z` (`GenRMT.norm_cformC_sub_cformC_le`, C3)
moves the point from the moving `z d / n` to the fixed `z'` deterministically (the bound does
not mention the matrix, so it holds for every outcome, not just in probability), the ratio
`d / n` is bounded by `Filter.Tendsto.bddAbove_range`, and `MP.c_mul_mTildeC_eq_mC_inv`
(`RMT/General/MPtilde.lean`) identifies the rescaled limit with `MP.mTildeC c z`. A three-term
split (moving the point, the random error at the fixed point, moving the coefficient) and
`TendstoInProb.add`/`const_mul`/`of_le` close L2.
-/

open Filter Topology MeasureTheory ProbabilityTheory Set
open scoped Matrix ENNReal NNReal

namespace StackedSVD
namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- Shared core of L1 and L2's random part: at a general noise matrix `Y` of shape `pN × dN`
with i.i.d. entries of law `ν`, along a regime `pN / dN → c'`, the quadratic form of the
resolvent at a sequence of unit vectors `u` converges to `MP.mC c' z` in probability. See the
route above the imports. -/
private theorem qform_tendsto_aux {pN dN : ℕ → ℕ} {c' : ℝ} (hc' : 0 < c')
    (hd : Tendsto dN atTop atTop) (hp : ∀ N, 0 < pN N) (hq : ∀ N, 0 < dN N)
    (hcN : Tendsto (fun N => (pN N : ℝ) / dN N) atTop (𝓝 c'))
    {ν : Measure ℝ} (hν : NoiseLaw ν)
    (Y : ∀ N, Ω N → Matrix (Fin (pN N)) (Fin (dN N)) ℝ)
    (hY : ∀ N, HasLaw (Y N) (noiseMatrix ν (pN N) (dN N)) (μ N))
    {z : ℂ} (hz : 0 < z.im) (u : ∀ N, Fin (dN N) → ℝ) (hu : ∀ N, u N ⬝ᵥ u N = 1) :
    TendstoInProb μ (fun N ω => ‖R4C.qformC (GenRMT.gram (Y N ω)) z (u N) - MP.mC c' z‖) 0 := by
  have hK : Tendsto (fun N => GenRMT.isoRate ν z (pN N) (dN N)) atTop (𝓝 0) :=
    GenRMT.tendsto_isoRate ν hz hd hcN
  have hbnd : ∀ ε : ℝ, 0 < ε → ∀ N,
      μ N {ω | ε ≤ |‖R4C.qformC (GenRMT.gram (Y N ω)) z (u N)
          - ((u N ⬝ᵥ u N : ℝ) : ℂ) * R4C.stieltjesC (GenRMT.gram (Y N ω)) z‖ - 0|}
        ≤ ENNReal.ofReal (GenRMT.isoRate ν z (pN N) (dN N) / ε ^ 2) := by
    intro ε hε N
    have hu1 : u N ⬝ᵥ u N ≤ 1 := (hu N).le
    have hset : MeasurableSet {Y' : Matrix (Fin (pN N)) (Fin (dN N)) ℝ | ε ≤
        ‖R4C.qformC (GenRMT.gram Y') z (u N)
          - ((u N ⬝ᵥ u N : ℝ) : ℂ) * R4C.stieltjesC (GenRMT.gram Y') z‖} :=
      measurableSet_le measurable_const (GenRMT.measurable_isoF z (u N)).norm
    have hbase := GenRMT.measure_qformC_sub_ge_le hν hz (hp N) (hq N) (u N) hu1 hε
    have heq := (hY N).measure_eq hset
    have hseteq : {ω | ε ≤ |‖R4C.qformC (GenRMT.gram (Y N ω)) z (u N)
          - ((u N ⬝ᵥ u N : ℝ) : ℂ) * R4C.stieltjesC (GenRMT.gram (Y N ω)) z‖ - 0|}
        = {ω | ε ≤ ‖R4C.qformC (GenRMT.gram (Y N ω)) z (u N)
            - ((u N ⬝ᵥ u N : ℝ) : ℂ) * R4C.stieltjesC (GenRMT.gram (Y N ω)) z‖} := by
      ext ω
      simp only [Set.mem_ofPred_eq, sub_zero, abs_norm]
    rw [hseteq, heq]
    exact hbase
  have hstep1raw := Cheb.tendstoInProb_of_meas_le hK hbnd
  have hcoef : ∀ N, ((u N ⬝ᵥ u N : ℝ) : ℂ) = 1 := fun N => by rw [hu N]; norm_num
  have hstep1 : TendstoInProb μ (fun N ω => ‖R4C.qformC (GenRMT.gram (Y N ω)) z (u N)
      - R4C.stieltjesC (GenRMT.gram (Y N ω)) z‖) 0 := by
    simpa only [hcoef, one_mul, sub_zero] using hstep1raw
  have hstep2 : TendstoInProb μ (fun N ω => ‖R4C.stieltjesC (GenRMT.gram (Y N ω)) z
      - MP.mC c' z‖) 0 :=
    GenRMT.tendstoInProb_stieltjesC_general hν hc' hz hd hp hcN Y hY
  have hsum := hstep1.add hstep2
  rw [add_zero] at hsum
  refine TendstoInProb.of_le (fun N => Filter.Eventually.of_forall fun ω => ?_) hsum
  rw [sub_zero, abs_of_nonneg (norm_nonneg _)]
  have hid : R4C.qformC (GenRMT.gram (Y N ω)) z (u N) - MP.mC c' z
      = (R4C.qformC (GenRMT.gram (Y N ω)) z (u N) - R4C.stieltjesC (GenRMT.gram (Y N ω)) z)
        + (R4C.stieltjesC (GenRMT.gram (Y N ω)) z - MP.mC c' z) := by ring
  rw [hid]
  exact norm_add_le _ _

/-! F37 (2026-09-09): `im_mul_natCast_div_natCast_pos'` and `tendstoInProb_of_tendsto_zero'`
used to be repeated here (`private`). `GenRMT.im_mul_natCast_div_natCast_pos`
(`FormsBridge.lean:389`, now public) and `SpikedModel.tendstoInProb_of_tendsto_zero`
(`R3minus.lean:398`, now public) are the canonical copies, used directly. -/

/-- **L1.** `vᵀ (gram Z - z)⁻¹ v → MP.mC c z` in probability: the isotropic quadratic form of
the resolvent at the deterministic direction `v` converges to the Marchenko-Pastur transform.
See the route above the imports. -/
theorem qformC_gram_v_tendsto [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hG : m.GeneralNoise ν) (z : ℂ) (hz : 0 < z.im) :
    TendstoInProb μ (fun N ω => ‖R4C.qformC (GenRMT.gram (m.Z N ω)) z
      (WithLp.ofLp (m.v N)) - MP.mC c z‖) 0 :=
  qform_tendsto_aux hc hreg.2.1 m.hn m.hd hreg.2.2 hν m.Z hG hz
    (fun N => WithLp.ofLp (m.v N)) (fun N => dotProduct_v_self m N)

/-- **L2.** `uᵀ (gramC Z - z)⁻¹ u → MP.mTildeC c z` in probability: the isotropic quadratic
form of the companion resolvent at the deterministic direction `u` converges to the companion
Marchenko-Pastur transform. See the route above the imports. -/
theorem qformC_gramC_u_tendsto [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hG : m.GeneralNoise ν) (z : ℂ) (hz : 0 < z.im) :
    TendstoInProb μ (fun N ω => ‖R4C.qformC (GenRMT.gramC (m.Z N ω)) z
      (WithLp.ofLp (m.u N)) - MP.mTildeC c z‖) 0 := by
  have : IsProbabilityMeasure ν := hν.prob
  -- the fixed point `z' = z / c`
  have hz' : 0 < (z / (c : ℂ)).im := by rw [Complex.div_ofReal_im]; exact div_pos hz hc
  have hdn : Tendsto (fun N => (d N : ℝ) / n N) atTop (𝓝 c⁻¹) := by
    have h := hreg.2.2.inv₀ hc.ne'
    simpa [inv_div] using h
  have hG' : ∀ N, HasLaw (fun ω => (m.Z N ω)ᵀ) (noiseMatrix ν (d N) (n N)) (μ N) := fun N =>
    (measurePreserving_transpose_noiseMatrix ν (n N) (d N)).fun_comp_hasLaw (hG N)
  -- Step A: the random part at the fixed point, on the transposed model
  have hFL : TendstoInProb μ (fun N ω => ‖R4C.qformC (GenRMT.gram ((m.Z N ω)ᵀ)) (z / (c : ℂ))
      (WithLp.ofLp (m.u N)) - MP.mC c⁻¹ (z / (c : ℂ))‖) 0 :=
    qform_tendsto_aux (inv_pos.mpr hc) hreg.1 m.hd m.hn hdn hν (fun N ω => (m.Z N ω)ᵀ) hG' hz'
      (fun N => WithLp.ofLp (m.u N)) (fun N => m.dotProduct_u_self N)
  -- the deterministic duality identity
  have hident : ∀ N ω, R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))
      = ((d N : ℂ) / (n N : ℂ)) * R4C.qformC (GenRMT.gram ((m.Z N ω)ᵀ)) (z * (d N : ℂ) / n N)
          (WithLp.ofLp (m.u N)) := fun N ω =>
    GenRMT.cformC_gramC_eq (m.Z N ω) (m.hn N) (m.hd N) hz _ _
  -- the moving-point Lipschitz bound, deterministic in `ω`
  have hpim : ∀ N, 0 < (z * (d N : ℂ) / (n N : ℂ)).im := fun N =>
    -- `GenRMT.im_mul_natCast_div_natCast_pos hz hp hd : 0 < (z * (d : ℂ) / (p : ℂ)).im`, so the
    -- denominator argument (`n N` here) comes first and the numerator argument (`d N`) second,
    -- the reverse of the deleted local `im_mul_natCast_div_natCast_pos'`.
    GenRMT.im_mul_natCast_div_natCast_pos hz (m.hn N) (m.hd N)
  have hLipPt : ∀ N ω, ‖R4C.qformC (GenRMT.gram ((m.Z N ω)ᵀ)) (z * (d N : ℂ) / n N)
        (WithLp.ofLp (m.u N))
      - R4C.qformC (GenRMT.gram ((m.Z N ω)ᵀ)) (z / (c : ℂ)) (WithLp.ofLp (m.u N))‖
      ≤ ‖z * (d N : ℂ) / n N - z / (c : ℂ)‖
          / ((z * (d N : ℂ) / n N).im * (z / (c : ℂ)).im) := by
    intro N ω
    have hb := GenRMT.norm_cformC_sub_cformC_le
      (GenRMT.gram_isHermitian ((m.Z N ω)ᵀ)) (hpim N) hz'
      (WithLp.ofLp (m.u N)) (WithLp.ofLp (m.u N))
    rwa [m.dotProduct_u_self N, Real.sqrt_one, mul_one, one_mul] at hb
  -- the point-moving bound vanishes deterministically
  have hpointTend : Tendsto (fun N => z * (d N : ℂ) / n N) atTop (𝓝 (z / (c : ℂ))) := by
    have h1 : Tendsto (fun N => (((d N : ℝ) / n N : ℝ) : ℂ)) atTop (𝓝 ((c⁻¹ : ℝ) : ℂ)) :=
      (Complex.continuous_ofReal.tendsto _).comp hdn
    have h2 := h1.const_mul z
    have hfun : ∀ N, z * (((d N : ℝ) / n N : ℝ) : ℂ) = z * (d N : ℂ) / n N := by
      intro N; push_cast; ring
    have hlim : z * ((c⁻¹ : ℝ) : ℂ) = z / (c : ℂ) := by
      rw [Complex.ofReal_inv, ← div_eq_mul_inv]
    rw [hlim] at h2
    exact h2.congr hfun
  have hLipTend : Tendsto (fun N => ‖z * (d N : ℂ) / n N - z / (c : ℂ)‖
        / ((z * (d N : ℂ) / n N).im * (z / (c : ℂ)).im)) atTop (𝓝 0) := by
    have hnum : Tendsto (fun N => ‖z * (d N : ℂ) / n N - z / (c : ℂ)‖) atTop (𝓝 0) := by
      have h : Tendsto (fun N => z * (d N : ℂ) / n N - z / (c : ℂ)) atTop
          (𝓝 (z / (c : ℂ) - z / (c : ℂ))) := hpointTend.sub tendsto_const_nhds
      rw [sub_self] at h
      simpa using h.norm
    have hden : Tendsto (fun N => (z * (d N : ℂ) / n N).im * (z / (c : ℂ)).im) atTop
        (𝓝 ((z / (c : ℂ)).im * (z / (c : ℂ)).im)) :=
      ((Complex.continuous_im.tendsto _).comp hpointTend).mul tendsto_const_nhds
    have hdenpos : (0 : ℝ) < (z / (c : ℂ)).im * (z / (c : ℂ)).im := mul_pos hz' hz'
    have hdiv := hnum.div hden hdenpos.ne'
    rwa [zero_div] at hdiv
  have hratTend : Tendsto (fun N => (d N : ℝ) / n N * (‖z * (d N : ℂ) / n N - z / (c : ℂ)‖
        / ((z * (d N : ℂ) / n N).im * (z / (c : ℂ)).im))) atTop (𝓝 0) := by
    have h := hdn.mul hLipTend
    simpa using h
  have hratNN : ∀ N, (0 : ℝ) ≤ (d N : ℝ) / n N * (‖z * (d N : ℂ) / n N - z / (c : ℂ)‖
        / ((z * (d N : ℂ) / n N).im * (z / (c : ℂ)).im)) := fun N => by
    refine mul_nonneg (by positivity) (div_nonneg (norm_nonneg _) ?_)
    exact (mul_pos (hpim N) hz').le
  have hDelta1 : TendstoInProb μ (fun N (_ : Ω N) => (d N : ℝ) / n N
      * (‖z * (d N : ℂ) / n N - z / (c : ℂ)‖
        / ((z * (d N : ℂ) / n N).im * (z / (c : ℂ)).im))) 0 :=
    tendstoInProb_of_tendsto_zero hratNN hratTend
  -- the coefficient `(d/n : ℂ)` is bounded (a convergent real sequence is bounded)
  obtain ⟨M, hM⟩ := hdn.bddAbove_range
  have hMr : ∀ N, (d N : ℝ) / n N ≤ M := fun N => hM (Set.mem_range_self N)
  have hMnn : (0 : ℝ) ≤ M := le_trans (by positivity : (0:ℝ) ≤ (d 0 : ℝ) / n 0) (hMr 0)
  have hDelta2 : TendstoInProb μ (fun N ω => M * ‖R4C.qformC
      (GenRMT.gram ((m.Z N ω)ᵀ)) (z / (c : ℂ)) (WithLp.ofLp (m.u N))
        - MP.mC c⁻¹ (z / (c : ℂ))‖) 0 := by
    have h := TendstoInProb.const_mul M hFL
    simpa using h
  -- the coefficient error `d/n - c⁻¹` vanishes deterministically
  have hrcTend : Tendsto (fun N => |(d N : ℝ) / n N - c⁻¹|) atTop (𝓝 0) := by
    have h := (hdn.sub_const c⁻¹).abs
    simpa using h
  have hDelta3Tend : Tendsto (fun N => |(d N : ℝ) / n N - c⁻¹|
      * ‖MP.mC c⁻¹ (z / (c : ℂ))‖) atTop (𝓝 0) := by
    have h := hrcTend.mul_const (‖MP.mC c⁻¹ (z / (c : ℂ))‖)
    simpa using h
  have hDelta3NN : ∀ N, (0 : ℝ) ≤ |(d N : ℝ) / n N - c⁻¹| * ‖MP.mC c⁻¹ (z / (c : ℂ))‖ :=
    fun N => by positivity
  have hDelta3 : TendstoInProb μ (fun N (_ : Ω N) => |(d N : ℝ) / n N - c⁻¹|
      * ‖MP.mC c⁻¹ (z / (c : ℂ))‖) 0 :=
    tendstoInProb_of_tendsto_zero hDelta3NN hDelta3Tend
  have hsum := (hDelta1.add hDelta2).add hDelta3
  rw [add_zero, add_zero] at hsum
  refine TendstoInProb.of_le (fun N => Filter.Eventually.of_forall fun ω => ?_) hsum
  rw [sub_zero, abs_of_nonneg (norm_nonneg _)]
  -- assemble: rewrite `mTildeC c z` and the duality identity, then split into three pieces
  have hcC : (c : ℂ) ≠ 0 := by exact_mod_cast hc.ne'
  have hmt : MP.mTildeC c z = (c : ℂ)⁻¹ * MP.mC c⁻¹ (z / (c : ℂ)) := by
    have hmul := MP.c_mul_mTildeC_eq_mC_inv hc hz
    rw [← hmul, ← mul_assoc, inv_mul_cancel₀ hcC, one_mul]
  have hccast : (c : ℂ)⁻¹ = ((c⁻¹ : ℝ) : ℂ) := (Complex.ofReal_inv c).symm
  rw [hident N ω, hmt, hccast]
  have hsplit : ((d N : ℂ) / (n N : ℂ))
        * R4C.qformC (GenRMT.gram ((m.Z N ω)ᵀ)) (z * (d N : ℂ) / n N) (WithLp.ofLp (m.u N))
      - ((c⁻¹ : ℝ) : ℂ) * MP.mC c⁻¹ (z / (c : ℂ))
      = ((d N : ℂ) / (n N : ℂ)) * (R4C.qformC (GenRMT.gram ((m.Z N ω)ᵀ))
            (z * (d N : ℂ) / n N) (WithLp.ofLp (m.u N))
          - R4C.qformC (GenRMT.gram ((m.Z N ω)ᵀ)) (z / (c : ℂ)) (WithLp.ofLp (m.u N)))
        + ((d N : ℂ) / (n N : ℂ)) * (R4C.qformC (GenRMT.gram ((m.Z N ω)ᵀ)) (z / (c : ℂ))
              (WithLp.ofLp (m.u N))
            - MP.mC c⁻¹ (z / (c : ℂ)))
        + (((d N : ℂ) / (n N : ℂ)) - ((c⁻¹ : ℝ) : ℂ)) * MP.mC c⁻¹ (z / (c : ℂ)) := by ring
  rw [hsplit]
  have hcoefeq : ((d N : ℂ) / (n N : ℂ)) = (((d N : ℝ) / n N : ℝ) : ℂ) := by push_cast; ring
  have hcoefnorm : ‖((d N : ℂ) / (n N : ℂ))‖ = (d N : ℝ) / n N := by
    rw [hcoefeq, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (by positivity)]
  have hcoefsubnorm : ‖((d N : ℂ) / (n N : ℂ)) - ((c⁻¹ : ℝ) : ℂ)‖
      = |(d N : ℝ) / n N - c⁻¹| := by
    rw [hcoefeq, ← Complex.ofReal_sub, Complex.norm_real, Real.norm_eq_abs]
  refine (norm_add_le _ _).trans (add_le_add (norm_add_le _ _) le_rfl) |>.trans ?_
  rw [norm_mul, norm_mul, norm_mul, hcoefnorm, hcoefsubnorm]
  refine add_le_add (add_le_add ?_ ?_) le_rfl
  · exact mul_le_mul_of_nonneg_left (hLipPt N ω) (by positivity)
  · exact mul_le_mul_of_nonneg_right (hMr N) (norm_nonneg _)

end SpikedModel
end StackedSVD
