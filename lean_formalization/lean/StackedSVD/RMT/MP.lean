/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Defs

/-!
# Item MP: the Marchenko-Pastur Stieltjes transform of `W₀`

Specification: `notes/archive/rmt_MP.md` (status `user OK` for proofs, 2026-08-29).
Numeric confirmation of every identity: `notes/archive/audit_numeric_2026-08-29.md` section C.

`m c` is a `def` given by the explicit formula, not an integral against `μ_c`
(modeling choice 1 of the note). Below the edge the formula returns junk, so every
lemma carries `bulkEdge c < z` or `bulkEdge c ≤ z` (modeling choice 2).

Consumers: R1, R3⁻, R5, R6, T.
-/

open Filter Topology

namespace StackedSVD
namespace MP

/-- Lower bulk edge `(1 - √c)²`. `bulkEdge c = (1 + √c)²` is in `Defs.lean`. -/
noncomputable def bulkEdgeLo (c : ℝ) : ℝ := (1 - Real.sqrt c) ^ 2

/-- Stieltjes transform of the MP law of `W₀`, real axis. Junk below the edge. -/
noncomputable def m (c z : ℝ) : ℝ :=
  (-(z + 1 - c) + Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c))) / (2 * z)

/-- Its derivative, in closed form. -/
noncomputable def mDeriv (c z : ℝ) : ℝ :=
  -(m c z * (m c z + 1)) / Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c))

section RealAxis

variable {c z θ : ℝ}

/-! ### Edge arithmetic -/

theorem bulkEdge_pos (_hc : 0 ≤ c) : 0 < bulkEdge c := by
  unfold bulkEdge
  have := Real.sqrt_nonneg c
  positivity

theorem bulkEdgeLo_lt_bulkEdge (hc : 0 < c) : bulkEdgeLo c < bulkEdge c := by
  unfold bulkEdge bulkEdgeLo
  have hs : 0 < Real.sqrt c := Real.sqrt_pos.mpr hc
  nlinarith

/-- MP2. -/
theorem discr_eq (hc : 0 ≤ c) :
    (z + 1 - c) ^ 2 - 4 * z = (z - bulkEdge c) * (z - bulkEdgeLo c) := by
  obtain ⟨s, hs0, rfl⟩ : ∃ s : ℝ, 0 ≤ s ∧ c = s ^ 2 :=
    ⟨Real.sqrt c, Real.sqrt_nonneg c, (Real.sq_sqrt hc).symm⟩
  unfold bulkEdge bulkEdgeLo
  rw [Real.sqrt_sq hs0]
  ring

theorem disc_nonneg (hc : 0 < c) (hz : bulkEdge c ≤ z) :
    0 ≤ (z - bulkEdge c) * (z - bulkEdgeLo c) := by
  have := bulkEdgeLo_lt_bulkEdge hc
  exact mul_nonneg (by linarith) (by linarith)

theorem disc_pos (hc : 0 < c) (hz : bulkEdge c < z) :
    0 < (z - bulkEdge c) * (z - bulkEdgeLo c) := by
  have := bulkEdgeLo_lt_bulkEdge hc
  exact mul_pos (by linarith) (by linarith)

theorem sqrt_disc_pos (hc : 0 < c) (hz : bulkEdge c < z) :
    0 < Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c)) :=
  Real.sqrt_pos.mpr (disc_pos hc hz)

theorem sq_sqrt_disc (hc : 0 < c) (hz : bulkEdge c ≤ z) :
    Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c)) ^ 2 = (z + 1 - c) ^ 2 - 4 * z := by
  rw [Real.sq_sqrt (disc_nonneg hc hz), ← discr_eq hc.le]

/-- Above the edge, `z + 1 - c ≥ 2 + 2√c > 0`. -/
theorem linear_pos (hc : 0 < c) (hz : bulkEdge c ≤ z) : 0 < z + 1 - c := by
  obtain ⟨s, hs0, rfl⟩ : ∃ s : ℝ, 0 < s ∧ c = s ^ 2 :=
    ⟨Real.sqrt c, Real.sqrt_pos.mpr hc, (Real.sq_sqrt hc.le).symm⟩
  unfold bulkEdge at hz
  rw [Real.sqrt_sq hs0.le] at hz
  nlinarith

/-- `bulkEdge c = 1 + 2√c + c`, the linear form used to avoid `√c ^ 2` bookkeeping. -/
theorem bulkEdge_eq (hc : 0 ≤ c) : bulkEdge c = 1 + 2 * Real.sqrt c + c := by
  unfold bulkEdge
  have hs2 : Real.sqrt c ^ 2 = c := Real.sq_sqrt hc
  linear_combination hs2

/-! ### The formula and the quadratic -/

/-- MP3.4. -/
theorem two_mul_mul_m_add (hc : 0 < c) (hz : bulkEdge c < z) :
    2 * z * m c z + z + 1 - c
      = Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c)) := by
  have hz0 : (0 : ℝ) < z := lt_of_lt_of_le (bulkEdge_pos hc.le) hz.le
  have hzne : z ≠ 0 := hz0.ne'
  have h : 2 * z * ((-(z + 1 - c) + Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c)))
      / (2 * z)) = -(z + 1 - c) + Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c)) := by
    field_simp
  unfold m
  rw [h]
  ring

/-- Rationalized form of `m`. Used for the sign, the lower bound and the limit at `∞`. -/
theorem m_eq_neg_two_div (hc : 0 < c) (hz : bulkEdge c < z) :
    m c z = -2 / (Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c)) + (z + 1 - c)) := by
  have hz0 : (0 : ℝ) < z := lt_of_lt_of_le (bulkEdge_pos hc.le) hz.le
  have hS := sqrt_disc_pos hc hz
  have hB := linear_pos hc hz.le
  have hS2 := sq_sqrt_disc hc hz.le
  unfold m
  rw [div_eq_div_iff (by linarith : (2 : ℝ) * z ≠ 0) (by linarith :
    Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c)) + (z + 1 - c) ≠ 0)]
  linear_combination hS2

/-- MP1. -/
theorem m_quadratic (hc : 0 < c) (hz : bulkEdge c < z) :
    z * m c z ^ 2 + (z + 1 - c) * m c z + 1 = 0 := by
  have hz0 : (0 : ℝ) < z := lt_of_lt_of_le (bulkEdge_pos hc.le) hz.le
  have hS2 := sq_sqrt_disc hc hz.le
  have h2 : 2 * z * m c z = Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c)) - (z + 1 - c) := by
    have := two_mul_mul_m_add hc hz; linarith
  have key : 4 * z * (z * m c z ^ 2 + (z + 1 - c) * m c z + 1) = 0 := by
    linear_combination
      (2 * z * m c z + (z + 1 - c) + Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c))) * h2 + hS2
  have h4 : (4 : ℝ) * z ≠ 0 := by positivity
  exact (mul_eq_zero.mp key).resolve_left h4

/-! ### Branch: `m` lands in `(-1/(1+√c), 0)` -/

theorem m_neg (hc : 0 < c) (hz : bulkEdge c < z) : m c z < 0 := by
  have hS := sqrt_disc_pos hc hz
  have hB := linear_pos hc hz.le
  rw [m_eq_neg_two_div hc hz]
  exact div_neg_of_neg_of_pos (by norm_num) (by linarith)

/-- The key positivity `0 < 1 + (1 + √c) * m c z`, equivalent to the lower bound of MP3.1. -/
theorem one_add_mul_m_pos (hc : 0 < c) (hz : bulkEdge c < z) :
    0 < 1 + (1 + Real.sqrt c) * m c z := by
  have hs2 : Real.sqrt c ^ 2 = c := Real.sq_sqrt hc.le
  have hS := sqrt_disc_pos hc hz
  have hB := linear_pos hc hz.le
  have hSBne : Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c)) + (z + 1 - c) ≠ 0 := by
    intro h; linarith [h ▸ (le_refl (0 : ℝ))]
  have hprod : (1 + (1 + Real.sqrt c) * m c z)
      * (Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c)) + (z + 1 - c))
      = Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c)) + (z + 1 - c)
          - 2 * (1 + Real.sqrt c) := by
    rw [m_eq_neg_two_div hc hz]
    field_simp
    ring
  have hshift : z + 1 - c - 2 - 2 * Real.sqrt c = z - bulkEdge c := by
    rw [bulkEdge_eq hc.le]; ring
  have hgap : 0 < Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c)) + (z + 1 - c)
      - 2 * (1 + Real.sqrt c) := by
    have : Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c)) + (z + 1 - c)
        - 2 * (1 + Real.sqrt c)
        = Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c)) + (z - bulkEdge c) := by
      linarith [hshift]
    rw [this]; linarith
  nlinarith [hprod, hgap, hB, hS]

/-- MP3.1. -/
theorem m_mem_Ioo (hc : 0 < c) (hz : bulkEdge c < z) :
    m c z ∈ Set.Ioo (-(1 / (1 + Real.sqrt c))) 0 := by
  have hs0 : 0 < Real.sqrt c := Real.sqrt_pos.mpr hc
  have h1s : (0 : ℝ) < 1 + Real.sqrt c := by linarith
  have h1sne : (1 : ℝ) + Real.sqrt c ≠ 0 := h1s.ne'
  refine ⟨?_, m_neg hc hz⟩
  have hpos := one_add_mul_m_pos hc hz
  have hdiff : m c z - -(1 / (1 + Real.sqrt c))
      = (1 + (1 + Real.sqrt c) * m c z) / (1 + Real.sqrt c) := by
    rw [eq_div_iff h1sne]
    field_simp
    ring
  have hd : 0 < m c z - -(1 / (1 + Real.sqrt c)) := by
    rw [hdiff]; exact div_pos hpos h1s
  linarith

/-- `-1 < m c z`, a corollary used for the sign of the derivative. -/
theorem neg_one_lt_m (hc : 0 < c) (hz : bulkEdge c < z) : -1 < m c z := by
  have hs0 : 0 < Real.sqrt c := Real.sqrt_pos.mpr hc
  have h1s : (0 : ℝ) < 1 + Real.sqrt c := by linarith
  have hlt : (1 : ℝ) / (1 + Real.sqrt c) < 1 := by
    rw [div_lt_one h1s]; linarith
  have := (m_mem_Ioo hc hz).1
  linarith

/-- MP3.3. -/
theorem m_bulkEdge (hc : 0 < c) : m c (bulkEdge c) = -(1 / (1 + Real.sqrt c)) := by
  obtain ⟨s, hs0, rfl⟩ : ∃ s : ℝ, 0 < s ∧ c = s ^ 2 :=
    ⟨Real.sqrt c, Real.sqrt_pos.mpr hc, (Real.sq_sqrt hc.le).symm⟩
  unfold m bulkEdge bulkEdgeLo
  simp only [Real.sqrt_sq hs0.le, sub_self, zero_mul, Real.sqrt_zero, add_zero]
  have h1s : (1 : ℝ) + s ≠ 0 := by positivity
  field_simp
  ring

/-! ### Continuity and the derivative -/

theorem continuousOn_m (hc : 0 < c) : ContinuousOn (m c) (Set.Ici (bulkEdge c)) := by
  have hb := bulkEdge_pos hc.le
  unfold m
  apply ContinuousOn.div
  · fun_prop
  · fun_prop
  · intro x hx
    simp only [Set.mem_Ici] at hx
    have hx0 : (0 : ℝ) < x := lt_of_lt_of_le hb hx
    intro h
    linarith [h ▸ (le_refl (0 : ℝ))]

private theorem deriv_algebra {z c S : ℝ} (hz0 : 0 < z) (hS0 : 0 < S)
    (hS2 : S ^ 2 = (z + 1 - c) ^ 2 - 4 * z) :
    -((-(z + 1 - c) + S) / (2 * z) * ((-(z + 1 - c) + S) / (2 * z) + 1)) / S
      = ((-1 + (2 * z - 2 - 2 * c) / (2 * S)) * (2 * z)
          - (-1 * z + (c - 1) + S) * 2) / (2 * z) ^ 2 := by
  have hzne : z ≠ 0 := hz0.ne'
  have hSne : S ≠ 0 := hS0.ne'
  have hden : (0 : ℝ) < 4 * z ^ 2 * S :=
    mul_pos (mul_pos (by norm_num) (pow_pos hz0 2)) hS0
  have hA : -((-(z + 1 - c) + S) / (2 * z) * ((-(z + 1 - c) + S) / (2 * z) + 1)) / S
      = (-((S - (z + 1 - c)) * (S - (z + 1 - c) + 2 * z))) / (4 * z ^ 2 * S) := by
    field_simp
    ring
  have hB : ((-1 + (2 * z - 2 - 2 * c) / (2 * S)) * (2 * z)
          - (-1 * z + (c - 1) + S) * 2) / (2 * z) ^ 2
      = (-2 * z * S + 2 * z ^ 2 - 2 * z - 2 * c * z - 2 * S ^ 2 + 2 * (z + 1 - c) * S)
          / (4 * z ^ 2 * S) := by
    field_simp
    ring
  rw [hA, hB, div_eq_div_iff hden.ne' hden.ne']
  linear_combination (4 * z ^ 2 * S) * hS2

theorem hasDerivAt_m (hc : 0 < c) (hz : bulkEdge c < z) :
    HasDerivAt (m c) (mDeriv c z) z := by
  have hz0 : (0 : ℝ) < z := lt_of_lt_of_le (bulkEdge_pos hc.le) hz.le
  have hzne : z ≠ 0 := hz0.ne'
  have hD := disc_pos hc hz
  have hS0 := sqrt_disc_pos hc hz
  have hS2 := sq_sqrt_disc hc hz.le
  have hs2 : Real.sqrt c ^ 2 = c := Real.sq_sqrt hc.le
  have h1 : HasDerivAt (fun w : ℝ => (w - bulkEdge c) * (w - bulkEdgeLo c))
      (2 * z - 2 - 2 * c) z := by
    have h := ((hasDerivAt_id' (x := z)).sub_const (bulkEdge c)).mul
      ((hasDerivAt_id' (x := z)).sub_const (bulkEdgeLo c))
    have e : (2 : ℝ) * z - 2 - 2 * c
        = 1 * (z - bulkEdgeLo c) + (z - bulkEdge c) * 1 := by
      unfold bulkEdge bulkEdgeLo
      linear_combination (2 : ℝ) * hs2
    rw [e]
    exact h
  have h2 := h1.sqrt hD.ne'
  have hlin : HasDerivAt (fun w : ℝ => -1 * w + (c - 1)) (-1) z := by
    have h : HasDerivAt (fun w : ℝ => -1 * w + (c - 1)) (-1 * 1) z :=
      ((hasDerivAt_id' (x := z)).const_mul (-1 : ℝ)).add_const (c - 1)
    rw [mul_one] at h
    exact h
  have h4 := hlin.add h2
  have h5 : HasDerivAt (fun w : ℝ => 2 * w) 2 z := by
    have h : HasDerivAt (fun w : ℝ => 2 * w) (2 * 1) z :=
      (hasDerivAt_id' (x := z)).const_mul (2 : ℝ)
    rw [mul_one] at h
    exact h
  have h2zne : (2 : ℝ) * z ≠ 0 := mul_ne_zero (by norm_num) hzne
  have h6 := h4.div h5 h2zne
  have hfun : m c = fun w : ℝ =>
      (-1 * w + (c - 1) + Real.sqrt ((w - bulkEdge c) * (w - bulkEdgeLo c))) / (2 * w) := by
    funext w
    unfold m
    ring
  have hval : mDeriv c z
      = ((-1 + (2 * z - 2 - 2 * c)
            / (2 * Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c)))) * (2 * z)
          - (-1 * z + (c - 1) + Real.sqrt ((z - bulkEdge c) * (z - bulkEdgeLo c))) * 2)
        / (2 * z) ^ 2 := by
    unfold mDeriv m
    exact deriv_algebra hz0 hS0 hS2
  rw [hfun, hval]
  exact h6

theorem mDeriv_pos (hc : 0 < c) (hz : bulkEdge c < z) : 0 < mDeriv c z := by
  have hS0 := sqrt_disc_pos hc hz
  have h1 := m_neg hc hz
  have h2 := neg_one_lt_m hc hz
  unfold mDeriv
  refine div_pos ?_ hS0
  nlinarith

/-- MP3.2. -/
theorem m_strictMonoOn (hc : 0 < c) : StrictMonoOn (m c) (Set.Ici (bulkEdge c)) := by
  apply strictMonoOn_of_deriv_pos (convex_Ici _) (continuousOn_m hc)
  intro x hx
  rw [interior_Ici] at hx
  simp only [Set.mem_Ioi] at hx
  rw [(hasDerivAt_m hc hx).deriv]
  exact mDeriv_pos hc hx

theorem m_tendsto_zero (hc : 0 < c) : Tendsto (m c) atTop (𝓝 0) := by
  have hEq : (m c) =ᶠ[atTop] fun x : ℝ =>
      -2 / (Real.sqrt ((x - bulkEdge c) * (x - bulkEdgeLo c)) + (x + 1 - c)) := by
    filter_upwards [eventually_gt_atTop (bulkEdge c)] with x hx using m_eq_neg_two_div hc hx
  have hB : Tendsto (fun x : ℝ => x + 1 - c) atTop atTop := by
    rw [tendsto_atTop_atTop]
    exact fun b => ⟨b + c, fun a ha => by linarith⟩
  have hden : Tendsto
      (fun x : ℝ => Real.sqrt ((x - bulkEdge c) * (x - bulkEdgeLo c)) + (x + 1 - c))
      atTop atTop := by
    refine tendsto_atTop_mono (fun x => ?_) hB
    have := Real.sqrt_nonneg ((x - bulkEdge c) * (x - bulkEdgeLo c))
    linarith
  have hinv := hden.inv_tendsto_atTop
  have hmain : Tendsto (fun x : ℝ =>
      -2 / (Real.sqrt ((x - bulkEdge c) * (x - bulkEdgeLo c)) + (x + 1 - c)))
      atTop (𝓝 0) := by
    have := hinv.const_mul (-2 : ℝ)
    simpa [Pi.inv_apply, div_eq_mul_inv] using this
  exact hmain.congr' hEq.symm

theorem continuousOn_mDeriv (hc : 0 < c) :
    ContinuousOn (mDeriv c) (Set.Ioi (bulkEdge c)) := by
  have hm : ContinuousOn (m c) (Set.Ioi (bulkEdge c)) :=
    (continuousOn_m hc).mono Set.Ioi_subset_Ici_self
  unfold mDeriv
  apply ContinuousOn.div
  · exact (hm.mul (hm.add continuousOn_const)).neg
  · fun_prop
  · intro x hx
    simp only [Set.mem_Ioi] at hx
    exact (sqrt_disc_pos hc hx).ne'

/-- MP3.4, second half. Needed by R3⁻. -/
theorem mDeriv_tendsto_atTop (hc : 0 < c) :
    Tendsto (mDeriv c) (𝓝[>] (bulkEdge c)) atTop := by
  have hs0 : 0 < Real.sqrt c := Real.sqrt_pos.mpr hc
  have hmt : Tendsto (m c) (𝓝[>] (bulkEdge c)) (𝓝 (m c (bulkEdge c))) := by
    have h1 : ContinuousWithinAt (m c) (Set.Ici (bulkEdge c)) (bulkEdge c) :=
      (continuousOn_m hc) _ (Set.mem_Ici.mpr le_rfl)
    exact h1.tendsto.mono_left (nhdsWithin_mono _ Set.Ioi_subset_Ici_self)
  have hnum : Tendsto (fun x => -(m c x * (m c x + 1))) (𝓝[>] (bulkEdge c))
      (𝓝 (-(m c (bulkEdge c) * (m c (bulkEdge c) + 1)))) :=
    (hmt.mul (hmt.add tendsto_const_nhds)).neg
  have hL : 0 < -(m c (bulkEdge c) * (m c (bulkEdge c) + 1)) := by
    rw [m_bulkEdge hc]
    have h1s : (0 : ℝ) < 1 + Real.sqrt c := by linarith
    have h2 : (0 : ℝ) < 1 / (1 + Real.sqrt c) := by positivity
    have h3 : (1 : ℝ) / (1 + Real.sqrt c) < 1 := by rw [div_lt_one h1s]; linarith
    nlinarith
  have hsq : Tendsto (fun x : ℝ => Real.sqrt ((x - bulkEdge c) * (x - bulkEdgeLo c)))
      (𝓝[>] (bulkEdge c)) (𝓝[>] 0) := by
    apply tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within
    · have hcont : Continuous (fun x : ℝ => (x - bulkEdge c) * (x - bulkEdgeLo c)) := by
        fun_prop
      have h0 : Tendsto (fun x : ℝ => (x - bulkEdge c) * (x - bulkEdgeLo c))
          (𝓝[>] (bulkEdge c)) (𝓝 0) := by
        have h : Tendsto (fun x : ℝ => (x - bulkEdge c) * (x - bulkEdgeLo c))
            (𝓝[>] (bulkEdge c))
            (𝓝 ((bulkEdge c - bulkEdge c) * (bulkEdge c - bulkEdgeLo c))) :=
          (hcont.tendsto (bulkEdge c)).mono_left nhdsWithin_le_nhds
        simpa using h
      simpa using h0.sqrt
    · filter_upwards [self_mem_nhdsWithin] with x hx
      exact sqrt_disc_pos hc hx
  have hinv := hsq.inv_tendsto_nhdsGT_zero
  have hfun : mDeriv c = fun x : ℝ => -(m c x * (m c x + 1))
      * (Real.sqrt ((x - bulkEdge c) * (x - bulkEdgeLo c)))⁻¹ := by
    funext x
    simp only [mDeriv, div_eq_mul_inv]
  rw [hfun]
  exact hnum.pos_mul_atTop hL hinv

/-! ### The outlier (MP4) -/

/-- MP4.2. -/
theorem rhoSq_eq (_hc : 0 < c) (hθ : 0 < θ) (h : c < θ ^ 4) :
    rhoSq θ c = (θ ^ 2 + 1) * (θ ^ 2 + c) / θ ^ 2 := by
  have hθne : θ ≠ 0 := hθ.ne'
  rw [rhoSq, if_pos h]
  field_simp
  ring

/-- MP4.3. -/
theorem rhoSq_sub_bulkEdge (hc : 0 < c) (hθ : 0 < θ) (h : c < θ ^ 4) :
    rhoSq θ c - bulkEdge c = (θ ^ 2 - Real.sqrt c) ^ 2 / θ ^ 2 := by
  have hθne : θ ≠ 0 := hθ.ne'
  rw [rhoSq_eq hc hθ h]
  obtain ⟨s, hs0, rfl⟩ : ∃ s : ℝ, 0 < s ∧ c = s ^ 2 :=
    ⟨Real.sqrt c, Real.sqrt_pos.mpr hc, (Real.sq_sqrt hc.le).symm⟩
  unfold bulkEdge
  rw [Real.sqrt_sq hs0.le]
  field_simp
  ring

theorem bulkEdge_lt_rhoSq (hc : 0 < c) (hθ : 0 < θ) (h : c < θ ^ 4) :
    bulkEdge c < rhoSq θ c := by
  have hne : θ ^ 2 - Real.sqrt c ≠ 0 := by
    intro hzero
    have h1 : Real.sqrt c = θ ^ 2 := by linarith
    have h2 : Real.sqrt c ^ 2 = c := Real.sq_sqrt hc.le
    rw [h1] at h2
    have h3 : c = θ ^ 4 := by rw [← h2]; ring
    linarith
  have hsq : 0 < (θ ^ 2 - Real.sqrt c) ^ 2 :=
    lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hne))
  have hpos : 0 < (θ ^ 2 - Real.sqrt c) ^ 2 / θ ^ 2 := div_pos hsq (pow_pos hθ 2)
  have := rhoSq_sub_bulkEdge hc hθ h
  linarith

/-- `√((ρ - b)(ρ - b')) = (θ⁴ - c)/θ²`. -/
theorem sqrt_disc_rhoSq (hc : 0 < c) (hθ : 0 < θ) (h : c < θ ^ 4) :
    Real.sqrt ((rhoSq θ c - bulkEdge c) * (rhoSq θ c - bulkEdgeLo c))
      = (θ ^ 4 - c) / θ ^ 2 := by
  have hθne : θ ≠ 0 := hθ.ne'
  have hprod : (rhoSq θ c - bulkEdge c) * (rhoSq θ c - bulkEdgeLo c)
      = ((θ ^ 4 - c) / θ ^ 2) ^ 2 := by
    rw [rhoSq_eq hc hθ h]
    obtain ⟨s, hs0, rfl⟩ : ∃ s : ℝ, 0 < s ∧ c = s ^ 2 :=
      ⟨Real.sqrt c, Real.sqrt_pos.mpr hc, (Real.sq_sqrt hc.le).symm⟩
    unfold bulkEdge bulkEdgeLo
    rw [Real.sqrt_sq hs0.le]
    field_simp
    ring
  rw [hprod, Real.sqrt_sq (div_nonneg (by linarith) (by positivity))]

/-- MP4.1, the value of `m` at the outlier. -/
theorem m_rhoSq (hc : 0 < c) (hθ : 0 < θ) (h : c < θ ^ 4) :
    m c (rhoSq θ c) = -(1 / (θ ^ 2 + 1)) := by
  have hθne : θ ≠ 0 := hθ.ne'
  have h1 : (θ : ℝ) ^ 2 + 1 ≠ 0 := by positivity
  have h2 : (θ : ℝ) ^ 2 + c ≠ 0 := by positivity
  unfold m
  rw [sqrt_disc_rhoSq hc hθ h, rhoSq_eq hc hθ h]
  field_simp
  ring

/-- MP4.4. -/
theorem mDeriv_rhoSq (hc : 0 < c) (hθ : 0 < θ) (h : c < θ ^ 4) :
    mDeriv c (rhoSq θ c) = θ ^ 4 / ((θ ^ 2 + 1) ^ 2 * (θ ^ 4 - c)) := by
  have hθne : θ ≠ 0 := hθ.ne'
  have h1 : (θ : ℝ) ^ 2 + 1 ≠ 0 := by positivity
  have h4 : (θ : ℝ) ^ 4 - c ≠ 0 := by intro hzero; linarith
  unfold mDeriv
  rw [m_rhoSq hc hθ h, sqrt_disc_rhoSq hc hθ h]
  field_simp
  ring

/-- MP4.5, the overlap identity. -/
theorem overlap_identity (hc : 0 < c) (hθ : 0 < θ) (h : c < θ ^ 4) :
    θ ^ 2 * m c (rhoSq θ c) ^ 2 / ((θ ^ 2 + 1) * mDeriv c (rhoSq θ c))
      = betaSq θ c := by
  have hθne : θ ≠ 0 := hθ.ne'
  have h1 : (θ : ℝ) ^ 2 + 1 ≠ 0 := by positivity
  have h4 : (θ : ℝ) ^ 4 - c ≠ 0 := by intro hzero; linarith
  have h5 : (θ : ℝ) ^ 4 + θ ^ 2 ≠ 0 := by positivity
  rw [m_rhoSq hc hθ h, mDeriv_rhoSq hc hθ h, betaSq, if_pos h]
  field_simp

/-- Branch lemma (MP4.1). The root `-1/(θ²+1)` of the secular equation lies in the
physical range of `m` exactly when `θ⁴ > c`. -/
theorem neg_inv_mem_Ioo_iff (hc : 0 < c) (hθ : 0 < θ) :
    -(1 / (θ ^ 2 + 1)) ∈ Set.Ioo (-(1 / (1 + Real.sqrt c))) 0 ↔ c < θ ^ 4 := by
  have hs0 : 0 < Real.sqrt c := Real.sqrt_pos.mpr hc
  have hs2 : Real.sqrt c ^ 2 = c := Real.sq_sqrt hc.le
  have h1s : (0 : ℝ) < 1 + Real.sqrt c := by linarith
  have ht : (0 : ℝ) < θ ^ 2 + 1 := by positivity
  have hneg : -(1 / (θ ^ 2 + 1)) < 0 := by
    have : (0 : ℝ) < 1 / (θ ^ 2 + 1) := by positivity
    linarith
  constructor
  · rintro ⟨hlt, -⟩
    have hst : Real.sqrt c < θ ^ 2 := by
      by_contra hcon
      have hcon' : θ ^ 2 ≤ Real.sqrt c := not_lt.mp hcon
      have hle : (1 : ℝ) / (1 + Real.sqrt c) ≤ 1 / (θ ^ 2 + 1) :=
        one_div_le_one_div_of_le ht (by linarith)
      linarith
    nlinarith
  · intro h
    refine ⟨?_, hneg⟩
    have hst : Real.sqrt c < θ ^ 2 := by
      rw [Real.sqrt_lt' (pow_pos hθ 2)]
      nlinarith
    have : (1 : ℝ) / (θ ^ 2 + 1) < 1 / (1 + Real.sqrt c) :=
      one_div_lt_one_div_of_lt h1s (by linarith)
    linarith

/-- MP4, uniqueness of the limiting secular root above the edge. -/
theorem secular_eq_zero_iff (hc : 0 < c) (hθ : 0 < θ) (h : c < θ ^ 4)
    (hz : bulkEdge c < z) : 1 + (θ ^ 2 + 1) * m c z = 0 ↔ z = rhoSq θ c := by
  have h1 : (θ : ℝ) ^ 2 + 1 ≠ 0 := by positivity
  constructor
  · intro heq
    have hmz : m c z = -(1 / (θ ^ 2 + 1)) := by
      field_simp
      linarith
    have hmρ : m c (rhoSq θ c) = -(1 / (θ ^ 2 + 1)) := m_rhoSq hc hθ h
    refine (m_strictMonoOn hc).injOn (Set.mem_Ici.mpr hz.le)
      (Set.mem_Ici.mpr (bulkEdge_lt_rhoSq hc hθ h).le) ?_
    rw [hmz, hmρ]
  · rintro rfl
    rw [m_rhoSq hc hθ h]
    field_simp
    norm_num

/-- MP5. This corrects the sign in `notes/archive/rmt_roadmap.md` row R6. -/
theorem secular_pos_of_subcritical (hc : 0 < c) (_hθ : 0 ≤ θ)
    (h : θ ^ 4 ≤ c) (hz : bulkEdge c < z) : 0 < 1 + (θ ^ 2 + 1) * m c z := by
  have hs0 : 0 < Real.sqrt c := Real.sqrt_pos.mpr hc
  have h1s : (0 : ℝ) < 1 + Real.sqrt c := by linarith
  have hst : θ ^ 2 ≤ Real.sqrt c := by
    rw [Real.le_sqrt (by positivity) hc.le]
    nlinarith
  have hlow := (m_mem_Ioo hc hz).1
  have hmul : (θ ^ 2 + 1) * (-(1 / (1 + Real.sqrt c))) < (θ ^ 2 + 1) * m c z :=
    mul_lt_mul_of_pos_left hlow (by positivity)
  have hfrac : (θ ^ 2 + 1) * (1 / (1 + Real.sqrt c)) ≤ 1 := by
    rw [mul_one_div, div_le_one h1s]
    linarith
  have hrw : (θ ^ 2 + 1) * (-(1 / (1 + Real.sqrt c)))
      = -((θ ^ 2 + 1) * (1 / (1 + Real.sqrt c))) := by ring
  rw [hrw] at hmul
  linarith

end RealAxis

/-! ### The complex group (MP6). Items R1 and T use this and nothing else.

Modeling choice 5 of the note: the statements are about the quadratic `quad`, with no
named complex branch of `m`. `vieta` is stated with a complex `z` (the note writes a real
`z`, which `existsUnique_root_im_pos` cannot use). -/

section ComplexPlane

/-- The MP quadratic with complex argument. -/
noncomputable def quad (c : ℝ) (z w : ℂ) : ℂ := z * w ^ 2 + (z + 1 - c) * w + 1

theorem quad_zero_ne_zero (c : ℝ) (z : ℂ) : quad c z 0 ≠ 0 := by
  simp [quad]

theorem root_ne_zero {c : ℝ} {z w : ℂ} (h : quad c z w = 0) : w ≠ 0 := by
  rintro rfl
  exact quad_zero_ne_zero c z h

/-- Vieta plus the substitution `P = 1 + u⁻¹`, `Q = 1 + v⁻¹`, which turns the quadratic
into `P Q = c`. -/
theorem vieta {c : ℝ} {z u v : ℂ} (hz : z ≠ 0) (hu : quad c z u = 0)
    (hv : quad c z v = 0) (huv : u ≠ v) :
    u * v = z⁻¹ ∧ (1 + u⁻¹) * (1 + v⁻¹) = (c : ℂ) := by
  have hune : u ≠ 0 := root_ne_zero hu
  have hvne : v ≠ 0 := root_ne_zero hv
  unfold quad at hu hv
  have hsum : z * (u + v) + (z + 1 - (c : ℂ)) = 0 := by
    have hfac : (u - v) * (z * (u + v) + (z + 1 - (c : ℂ))) = 0 := by
      linear_combination hu - hv
    rcases mul_eq_zero.mp hfac with h | h
    · exact absurd (sub_eq_zero.mp h) huv
    · exact h
  have hprod : z * (u * v) = 1 := by
    linear_combination u * hsum - hu
  refine ⟨?_, ?_⟩
  · field_simp
    linear_combination hprod
  · have hkey : (u + 1) * (v + 1) = (c : ℂ) * (u * v) :=
      mul_left_cancel₀ hz (by linear_combination (1 - (c : ℂ)) * hprod + hsum)
    field_simp
    linear_combination hkey

/-- Real core of the branch lemma. If `P Q = c > 0` with both imaginary parts `≤ 0`, both
imaginary parts are `0`. -/
private theorem im_zero_core {c pr pi qr qi : ℝ} (hc : 0 < c)
    (e1 : pr * qr - pi * qi = c) (e2 : pr * qi + pi * qr = 0)
    (hp : pi ≤ 0) (hq : qi ≤ 0) : pi = 0 ∧ qi = 0 := by
  have key1 : (pr ^ 2 + pi ^ 2) * qi = -c * pi := by linear_combination pr * e2 - pi * e1
  have key2 : (qr ^ 2 + qi ^ 2) * pi = -c * qi := by linear_combination qr * e2 - qi * e1
  have h1 : (pr ^ 2 + pi ^ 2) * qi ≤ 0 := by nlinarith [sq_nonneg pr, sq_nonneg pi]
  have h2 : (qr ^ 2 + qi ^ 2) * pi ≤ 0 := by nlinarith [sq_nonneg qr, sq_nonneg qi]
  have hp0 : 0 ≤ pi := by nlinarith [key1, h1]
  have hq0 : 0 ≤ qi := by nlinarith [key2, h2]
  exact ⟨le_antisymm hp hp0, le_antisymm hq hq0⟩

private theorem re_im_eqs {c : ℝ} {u v : ℂ} (h : (1 + u⁻¹) * (1 + v⁻¹) = (c : ℂ)) :
    (1 + u⁻¹).re * (1 + v⁻¹).re - (1 + u⁻¹).im * (1 + v⁻¹).im = c ∧
      (1 + u⁻¹).re * (1 + v⁻¹).im + (1 + u⁻¹).im * (1 + v⁻¹).re = 0 := by
  constructor
  · rw [← Complex.mul_re, h, Complex.ofReal_re]
  · rw [← Complex.mul_im, h, Complex.ofReal_im]

private theorem im_inv_add_one (u : ℂ) : (1 + u⁻¹).im = -u.im / Complex.normSq u := by
  simp [Complex.inv_im]

/-- Branch lemma, upper half plane. Two roots of the MP quadratic cannot both have
`Im ≥ 0` unless both are real. -/
theorem im_eq_zero_of_im_nonneg {c : ℝ} (hc : 0 < c) {u v : ℂ} (hu : u ≠ 0) (hv : v ≠ 0)
    (h : (1 + u⁻¹) * (1 + v⁻¹) = (c : ℂ)) (hu' : 0 ≤ u.im) (hv' : 0 ≤ v.im) :
    u.im = 0 ∧ v.im = 0 := by
  have hnu : 0 < Complex.normSq u := Complex.normSq_pos.mpr hu
  have hnv : 0 < Complex.normSq v := Complex.normSq_pos.mpr hv
  obtain ⟨e1, e2⟩ := re_im_eqs h
  have hp : (1 + u⁻¹).im ≤ 0 := by
    rw [im_inv_add_one]
    exact div_nonpos_of_nonpos_of_nonneg (by linarith) hnu.le
  have hq : (1 + v⁻¹).im ≤ 0 := by
    rw [im_inv_add_one]
    exact div_nonpos_of_nonpos_of_nonneg (by linarith) hnv.le
  obtain ⟨hpz, hqz⟩ := im_zero_core hc e1 e2 hp hq
  rw [im_inv_add_one] at hpz hqz
  constructor
  · rcases div_eq_zero_iff.mp hpz with h' | h'
    · linarith
    · exact absurd h' hnu.ne'
  · rcases div_eq_zero_iff.mp hqz with h' | h'
    · linarith
    · exact absurd h' hnv.ne'

/-- Branch lemma, lower half plane (the mirror of `im_eq_zero_of_im_nonneg`). -/
theorem im_eq_zero_of_im_nonpos {c : ℝ} (hc : 0 < c) {u v : ℂ} (hu : u ≠ 0) (hv : v ≠ 0)
    (h : (1 + u⁻¹) * (1 + v⁻¹) = (c : ℂ)) (hu' : u.im ≤ 0) (hv' : v.im ≤ 0) :
    u.im = 0 ∧ v.im = 0 := by
  have hnu : 0 < Complex.normSq u := Complex.normSq_pos.mpr hu
  have hnv : 0 < Complex.normSq v := Complex.normSq_pos.mpr hv
  obtain ⟨e1, e2⟩ := re_im_eqs h
  have e1' : (1 + u⁻¹).re * (1 + v⁻¹).re - (-(1 + u⁻¹).im) * (-(1 + v⁻¹).im) = c := by
    linear_combination e1
  have e2' : (1 + u⁻¹).re * (-(1 + v⁻¹).im) + (-(1 + u⁻¹).im) * (1 + v⁻¹).re = 0 := by
    linear_combination -e2
  have hp : -(1 + u⁻¹).im ≤ 0 := by
    rw [im_inv_add_one]
    have : 0 ≤ -u.im / Complex.normSq u := div_nonneg (by linarith) hnu.le
    linarith
  have hq : -(1 + v⁻¹).im ≤ 0 := by
    rw [im_inv_add_one]
    have : 0 ≤ -v.im / Complex.normSq v := div_nonneg (by linarith) hnv.le
    linarith
  obtain ⟨hpz, hqz⟩ := im_zero_core hc e1' e2' hp hq
  have hpz' : (1 + u⁻¹).im = 0 := by linarith
  have hqz' : (1 + v⁻¹).im = 0 := by linarith
  rw [im_inv_add_one] at hpz' hqz'
  constructor
  · rcases div_eq_zero_iff.mp hpz' with h' | h'
    · linarith
    · exact absurd h' hnu.ne'
  · rcases div_eq_zero_iff.mp hqz' with h' | h'
    · linarith
    · exact absurd h' hnv.ne'

/-- Either square root of the discriminant gives a root of `quad`. -/
theorem quad_root_of_sq {c : ℝ} {z e : ℂ} (hz : z ≠ 0)
    (he : e ^ 2 = (z + 1 - (c : ℂ)) ^ 2 - 4 * z) :
    quad c z ((-(z + 1 - (c : ℂ)) + e) / (2 * z)) = 0 := by
  have h2z : ((2 : ℂ) * z) ≠ 0 := mul_ne_zero (by norm_num) hz
  have key : ((2 : ℂ) * z) ^ 2 * quad c z ((-(z + 1 - (c : ℂ)) + e) / (2 * z))
      = z * (e ^ 2 - ((z + 1 - (c : ℂ)) ^ 2 - 4 * z)) := by
    unfold quad
    field_simp
    ring
  rw [he, sub_self, mul_zero] at key
  exact (mul_eq_zero.mp key).resolve_left (pow_ne_zero 2 h2z)

/-- MP6. This is the statement R1 uses to identify its limit. -/
theorem existsUnique_root_im_pos {c : ℝ} (hc : 0 < c) {z : ℂ} (hz : 0 < z.im) :
    ∃! w : ℂ, 0 < w.im ∧ quad c z w = 0 := by
  have hz0 : z ≠ 0 := by
    intro h
    rw [h] at hz
    simp at hz
  have huniq : ∀ x y : ℂ, (0 < x.im ∧ quad c z x = 0) → (0 < y.im ∧ quad c z y = 0) →
      x = y := by
    rintro x y ⟨hx1, hx2⟩ ⟨hy1, hy2⟩
    by_contra hne
    obtain ⟨-, hPQ⟩ := vieta hz0 hx2 hy2 hne
    obtain ⟨hxi, -⟩ :=
      im_eq_zero_of_im_nonneg hc (root_ne_zero hx2) (root_ne_zero hy2) hPQ hx1.le hy1.le
    linarith
  obtain ⟨d, hd⟩ : ∃ d : ℂ, d ^ 2 = (z + 1 - (c : ℂ)) ^ 2 - 4 * z :=
    IsAlgClosed.exists_pow_nat_eq _ (by norm_num)
  have hd' : (-d) ^ 2 = (z + 1 - (c : ℂ)) ^ 2 - 4 * z := by rw [neg_sq]; exact hd
  have hr1 := quad_root_of_sq (c := c) hz0 hd
  have hr2 := quad_root_of_sq (c := c) hz0 hd'
  set w1 : ℂ := (-(z + 1 - (c : ℂ)) + d) / (2 * z) with hw1
  set w2 : ℂ := (-(z + 1 - (c : ℂ)) + -d) / (2 * z) with hw2
  have hn1 : w1 ≠ 0 := root_ne_zero hr1
  have hn2 : w2 ≠ 0 := root_ne_zero hr2
  have hsum : w1 + w2 = -(z + 1 - (c : ℂ)) / z := by
    rw [hw1, hw2]
    field_simp
    ring
  have hprod : w1 * w2 = z⁻¹ := by
    have hkey : ((2 : ℂ) * z) ^ 2 * (w1 * w2) = ((2 : ℂ) * z) ^ 2 * z⁻¹ := by
      rw [hw1, hw2]
      field_simp
      linear_combination -hd
    exact mul_left_cancel₀ (pow_ne_zero 2 (mul_ne_zero (by norm_num) hz0)) hkey
  have hPQ : (1 + w1⁻¹) * (1 + w2⁻¹) = (c : ℂ) := by
    have hkey : (w1 + 1) * (w2 + 1) = (c : ℂ) * (w1 * w2) := by
      have e : (w1 + 1) * (w2 + 1) = w1 * w2 + (w1 + w2) + 1 := by ring
      rw [e, hsum, hprod]
      field_simp
      ring
    field_simp
    linear_combination hkey
  have hone : 0 < w1.im ∨ 0 < w2.im := by
    rcases lt_or_ge 0 w1.im with h | h1
    · exact Or.inl h
    rcases lt_or_ge 0 w2.im with h | h2
    · exact Or.inr h
    exfalso
    obtain ⟨hi1, hi2⟩ := im_eq_zero_of_im_nonpos hc hn1 hn2 hPQ h1 h2
    have hzi : (z⁻¹).im = 0 := by
      rw [← hprod, Complex.mul_im, hi1, hi2]
      ring
    rw [Complex.inv_im] at hzi
    have hnz : 0 < Complex.normSq z := Complex.normSq_pos.mpr hz0
    rcases div_eq_zero_iff.mp hzi with h' | h'
    · linarith
    · exact absurd h' hnz.ne'
  rcases hone with h | h
  · exact ⟨w1, ⟨h, hr1⟩, fun y hy => huniq y w1 hy ⟨h, hr1⟩⟩
  · exact ⟨w2, ⟨h, hr2⟩, fun y hy => huniq y w2 hy ⟨h, hr2⟩⟩

/-- Link to the real axis, for item T. -/
theorem quad_ofReal {c z : ℝ} (hc : 0 < c) (hz : bulkEdge c < z) :
    quad c (z : ℂ) ((m c z : ℝ) : ℂ) = 0 := by
  have h := m_quadratic hc hz
  have hcast : ((z * m c z ^ 2 + (z + 1 - c) * m c z + 1 : ℝ) : ℂ) = 0 := by
    rw [h]
    simp
  unfold quad
  push_cast at hcast
  linear_combination hcast

end ComplexPlane

end MP
end StackedSVD
