/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Deloc
import StackedSVD.RMT.General.DelocLimits
import StackedSVD.RMT.General.R3minus
import StackedSVD.RMT.General.Simplicity

/-!
# Item R6' at a general law: the subcritical `align` field

The assembly of `align_tendstoInProb_of_subcritical_general`, the frozen statement of the
plan (`notes/archive/prop_single_table_general.md`, section 2; held in
`RMT/General/Statements.lean` until its discharge, the file retired by F35).
The skeleton is the Gaussian proof `RMT/R6.lean:879` to `:990`; only the block that bounds the
component of `v` along `r = v - a q` changes. The Gaussian proof gets that block from rotation
invariance (`m.tendstoInProb_overlap_rdir`); at a general law the same bound comes from a
spectral window: the top eigenvalue of the Gram matrix sits inside `[E - η, E + η]` with
`E = bulkEdge c`, so `overlap X r` is at most the squared norm of the spectral projector of
that window, and that is at most `2 η Im (rᵀ G(E + i η) r)` (`Deloc.normSq_specProj_Icc_le`).
The complex form on the Gram matrix is a rank-one downdate of the one on `W₀`
(`Deloc.cformC_add_vecMulVec_eq`), whose limit `m(z)/(θ² + 1)` stays bounded as `η ↓ 0`
(`GenRMT.norm_mC_edge_le`), so `η` can be taken small.
-/

open Filter Topology MeasureTheory ProbabilityTheory
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

namespace GenRMT

/-! ## A bound on `mC` at the bulk edge, uniform in the height -/

/-- A crude upper bound for `‖MP.mC c (bulkEdge c + i η)‖`, valid for every `0 < η ≤ 1`.
Only its finiteness matters: it lets the window height `η` be chosen after it. -/
noncomputable def mCbound (c : ℝ) : ℝ :=
  (bulkEdge c + 2 + c) + (((bulkEdge c + 2 + c) ^ 2 + 4 * (bulkEdge c + 1)) + 1) / 2

theorem mCbound_nonneg {c : ℝ} (hc : 0 < c) : 0 ≤ mCbound c := by
  have h1 : (0 : ℝ) ≤ bulkEdge c := by
    rw [bulkEdge]; positivity
  rw [mCbound]
  nlinarith [sq_nonneg (bulkEdge c + 2 + c)]

/-- **The uniform bound.** `mC c z = (-(z + 1 - c) + √((z-b)(z-b')))/(2z)`; on the segment
`z = bulkEdge c + i η`, `0 < η ≤ 1`, the numerator is bounded by the polynomial `mCbound c`
and `‖2 z‖ ≥ 2`, because `bulkEdge c ≥ 1`. -/
theorem norm_mC_edge_le {c : ℝ} (hc : 0 < c) {η : ℝ} (hη0 : 0 < η) (hη1 : η ≤ 1) :
    ‖MP.mC c ((bulkEdge c : ℂ) + (η : ℂ) * Complex.I)‖ ≤ mCbound c := by
  set E : ℝ := bulkEdge c with hEdef
  set z : ℂ := (E : ℂ) + (η : ℂ) * Complex.I with hzdef
  have hE1 : 1 ≤ E := by
    rw [hEdef, bulkEdge]
    nlinarith [Real.sqrt_nonneg c]
  have hzre : z.re = E := by simp [hzdef]
  have hzim : z.im = η := by simp [hzdef]
  have hzimpos : 0 < z.im := by rw [hzim]; exact hη0
  have hzn1 : 1 ≤ ‖z‖ := by
    have h : |z.re| ≤ ‖z‖ := Complex.abs_re_le_norm z
    rw [hzre, abs_of_nonneg (by linarith : (0:ℝ) ≤ E)] at h
    linarith
  have hznle : ‖z‖ ≤ E + 1 := by
    have h : ‖z‖ ≤ ‖(E : ℂ)‖ + ‖(η : ℂ) * Complex.I‖ := by
      rw [hzdef]; exact norm_add_le _ _
    rw [Complex.norm_real, norm_mul, Complex.norm_I, Complex.norm_real, mul_one,
      Real.norm_eq_abs, Real.norm_eq_abs, abs_of_nonneg (by linarith : (0:ℝ) ≤ E),
      abs_of_nonneg hη0.le] at h
    linarith
  -- the linear part of the numerator
  have hT : ‖z + 1 - (c : ℂ)‖ ≤ E + 2 + c := by
    have h1 : ‖z + 1 - (c : ℂ)‖ ≤ ‖z + 1‖ + ‖(c : ℂ)‖ := norm_sub_le _ _
    have h2 : ‖z + 1‖ ≤ ‖z‖ + ‖(1 : ℂ)‖ := norm_add_le _ _
    rw [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hc.le] at h1
    rw [norm_one] at h2
    linarith
  -- the square root part
  have hsq : ‖MP.sqrtDisc c z‖ ^ 2 ≤ (E + 2 + c) ^ 2 + 4 * (E + 1) := by
    have he : ‖MP.sqrtDisc c z‖ ^ 2 = ‖(z + 1 - (c : ℂ)) ^ 2 - 4 * z‖ := by
      rw [← norm_pow, MP.sqrtDisc_sq hc.le hzimpos]
    have h1 : ‖(z + 1 - (c : ℂ)) ^ 2 - 4 * z‖ ≤ ‖(z + 1 - (c : ℂ)) ^ 2‖ + ‖(4 : ℂ) * z‖ :=
      norm_sub_le _ _
    have h2 : ‖(z + 1 - (c : ℂ)) ^ 2‖ = ‖z + 1 - (c : ℂ)‖ ^ 2 := norm_pow _ 2
    have h3 : ‖(4 : ℂ) * z‖ = 4 * ‖z‖ := by
      rw [norm_mul]; norm_num
    have h4 : ‖z + 1 - (c : ℂ)‖ ^ 2 ≤ (E + 2 + c) ^ 2 :=
      pow_le_pow_left₀ (norm_nonneg _) hT 2
    rw [he]
    linarith
  have hsd : ‖MP.sqrtDisc c z‖ ≤ (((E + 2 + c) ^ 2 + 4 * (E + 1)) + 1) / 2 := by
    nlinarith [sq_nonneg (‖MP.sqrtDisc c z‖ - 1), norm_nonneg (MP.sqrtDisc c z)]
  -- assemble
  have hnum : ‖-(z + 1 - (c : ℂ)) + MP.sqrtDisc c z‖ ≤ mCbound c := by
    have h : ‖-(z + 1 - (c : ℂ)) + MP.sqrtDisc c z‖
        ≤ ‖-(z + 1 - (c : ℂ))‖ + ‖MP.sqrtDisc c z‖ := norm_add_le _ _
    rw [norm_neg] at h
    rw [mCbound, ← hEdef]
    linarith
  have hmCdef : MP.mC c z = (-(z + 1 - (c : ℂ)) + MP.sqrtDisc c z) / (2 * z) := rfl
  have hden : (2 : ℝ) ≤ ‖(2 : ℂ) * z‖ := by
    rw [norm_mul]
    norm_num
    linarith
  rw [hmCdef, norm_div, div_le_iff₀ (by linarith : (0:ℝ) < ‖(2 : ℂ) * z‖)]
  nlinarith [mCbound_nonneg hc, hnum, hden]

end GenRMT

/-! ## The assembly -/

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- **Item R6' at a general law.** For `θ⁴ ≤ c` the overlap of the top right singular vector
of `X` with `v` vanishes in probability, which is the `align` field of `SingleTableLaw`
because `betaSq θ c = 0` there. Twin of `align_tendstoInProb_of_subcritical`
(`RMT/R6.lean:879`), whose rotation-invariance block is replaced by the spectral window
bound of `RMT/General/Deloc.lean`. -/
theorem align_tendstoInProb_of_subcritical_general' [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ}
    (hc : 0 < c) (m : SpikedModel μ n d) (H : m.ResolventFormsC c) (RL : m.ResolventLimits c)
    (hreg : m.Regime c) {ν : Measure ℝ} [SigmaFinite ν] (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.GeneralNoise ν) (hθ : m.θ ^ 4 ≤ c) :
    TendstoInProb μ (fun N ω => overlap (m.X N ω) (m.v N)) (betaSq m.θ c) := by
  have hbeta : betaSq m.θ c = 0 := by rw [betaSq, if_neg (not_lt.mpr hθ)]
  rw [hbeta]
  intro ε hε
  have hK : (0 : ℝ) < m.θ ^ 2 + 1 := by positivity
  -- the real point above the edge, for the `q` component
  obtain ⟨z₀, hz₀b, -, hz₀K⟩ :=
    GenRMT.R3minus.exists_mDeriv_gt hc (16 / ((m.θ ^ 2 + 1) ^ 2 * ε))
  have hmD : 0 < MP.mDeriv c z₀ := MP.mDeriv_pos hc hz₀b
  have hbound : 8 < ε / 2 * ((m.θ ^ 2 + 1) ^ 2 * MP.mDeriv c z₀) := by
    rw [div_lt_iff₀ (by positivity)] at hz₀K
    nlinarith
  -- the window height and the two complex tolerances, for the `r` component
  set Cb : ℝ := GenRMT.mCbound c with hCbdef
  have hCb0 : 0 ≤ Cb := GenRMT.mCbound_nonneg hc
  set S : ℝ := Real.sqrt (m.θ ^ 2 + 2) with hSdef
  have hS0 : 0 ≤ S := Real.sqrt_nonneg _
  set η : ℝ := min (min 1 (ε / (24 * (Cb + 1)))) (z₀ - bulkEdge c) with hηdef
  have hη0 : 0 < η := lt_min (lt_min one_pos (by positivity)) (by linarith)
  have hη1 : η ≤ 1 := le_trans (min_le_left _ _) (min_le_left _ _)
  have hηz : η ≤ z₀ - bulkEdge c := min_le_right _ _
  have hηC : 24 * (η * Cb) + 24 * η ≤ ε := by
    have h : η ≤ ε / (24 * (Cb + 1)) := le_trans (min_le_left _ _) (min_le_right _ _)
    have h2 : η * (24 * (Cb + 1)) ≤ ε := by
      rw [← le_div_iff₀ (by positivity)]
      exact h
    nlinarith [h2]
  set δ₁ : ℝ := ε / 48 with hδ₁def
  have hδ₁ : 0 < δ₁ := by rw [hδ₁def]; positivity
  set δ₂ : ℝ := ε / (24 * (S + 1)) with hδ₂def
  have hδ₂ : 0 < δ₂ := by rw [hδ₂def]; positivity
  have hδ₂e : 24 * (S * δ₂) + 24 * δ₂ = ε := by
    rw [hδ₂def]
    field_simp
  set z : ℂ := ((bulkEdge c : ℝ) : ℂ) + (η : ℂ) * Complex.I with hzdef
  have hzim : z.im = η := by simp [hzdef]
  have hzimpos : 0 < z.im := by rw [hzim]; exact hη0
  have hmCb : ‖MP.mC c z‖ ≤ Cb := by
    rw [hCbdef, hzdef]
    exact GenRMT.norm_mC_edge_le hc hη0 hη1
  -- the six bad families
  have ht2pos : (0 : ℝ) < min ((m.θ ^ 2 + 1) / 2) 1 := lt_min (by positivity) one_pos
  have hB2 := tendstoInProb_dotProduct_qvec_general m hν hG hreg.2.1
    (min ((m.θ ^ 2 + 1) / 2) 1) ht2pos
  have hB3 := GenRMT.Deloc.tendstoInProb_qform2_qvec m RL hz₀b
    ((m.θ ^ 2 + 1) * MP.mDeriv c z₀ / 2) (by positivity)
  have hrho : rhoSq m.θ c = bulkEdge c := by rw [rhoSq, if_neg (not_lt.mpr hθ)]
  have hB4 := lamMax_tendstoInProb_of_subcritical_general RL hc hreg hν hG hθ η hη0
  rw [hrho] at hB4
  have hB5 := GenRMT.Deloc.qformC_W0_rdir_tendsto m H hreg hν hG z hzimpos δ₁ hδ₁
  have hB6 := GenRMT.Deloc.cformC_W0_qvec_rdir_tendsto m H hreg hν hG z hzimpos δ₂ hδ₂
  have hB7 : Tendsto (fun N => μ N {ω : Ω N | ¬ TopSimple ((m.X N ω)ᵀ * m.X N ω)
      (isHermitian_transpose_mul_self (m.X N ω))}) atTop (𝓝 0) := by
    have hz : ∀ N, μ N {ω : Ω N | ¬ TopSimple ((m.X N ω)ᵀ * m.X N ω)
        (isHermitian_transpose_mul_self (m.X N ω))} = 0 := fun N =>
      ae_iff.mp (singleTableLaw_topSimple_of_general m hac hG N)
    simp only [hz]
    exact tendsto_const_nhds
  refine tendsto_measure_zero_of_subset ?_ (tendsto_measure_zero_union hB2
    (tendsto_measure_zero_union hB3 (tendsto_measure_zero_union hB4
      (tendsto_measure_zero_union hB5 (tendsto_measure_zero_union hB6 hB7)))))
  intro N ω hω
  by_contra hcon
  simp only [Set.mem_union, not_or] at hcon
  obtain ⟨g2, g3, g4, g5, g6, g7⟩ := hcon
  -- unpack the good event
  have hqq2 : |m.qvec N ω ⬝ᵥ m.qvec N ω - (m.θ ^ 2 + 1)| < min ((m.θ ^ 2 + 1) / 2) 1 := by
    have h : ¬ (min ((m.θ ^ 2 + 1) / 2) 1
        ≤ |m.qvec N ω ⬝ᵥ m.qvec N ω - (m.θ ^ 2 + 1)|) := g2
    exact not_le.mp h
  have hqqlo : (m.θ ^ 2 + 1) / 2 < m.qvec N ω ⬝ᵥ m.qvec N ω := by
    have h := (abs_lt.mp hqq2).1
    have h2 : min ((m.θ ^ 2 + 1) / 2) 1 ≤ (m.θ ^ 2 + 1) / 2 := min_le_left _ _
    linarith
  have hqqhi : m.qvec N ω ⬝ᵥ m.qvec N ω ≤ m.θ ^ 2 + 2 := by
    have h := (abs_lt.mp hqq2).2
    have h2 : min ((m.θ ^ 2 + 1) / 2) 1 ≤ 1 := min_le_right _ _
    linarith
  have hQ2 : (m.θ ^ 2 + 1) * MP.mDeriv c z₀ / 2
      < R4.qform2 (m.W0 N ω) z₀ (m.qvec N ω) := by
    have h0 : ¬ ((m.θ ^ 2 + 1) * MP.mDeriv c z₀ / 2
        ≤ |R4.qform2 (m.W0 N ω) z₀ (m.qvec N ω) - (m.θ ^ 2 + 1) * MP.mDeriv c z₀|) := g3
    have h := abs_lt.mp (not_le.mp h0)
    linarith [h.1]
  have hg4 : |gramLamMax (m.X N ω) - bulkEdge c| < η := by
    have h : ¬ (η ≤ |gramLamMax (m.X N ω) - bulkEdge c|) := g4
    exact not_le.mp h
  have hwin : |gramLamMax (m.X N ω) - bulkEdge c| ≤ η := hg4.le
  have hlam : gramLamMax (m.X N ω) ≤ z₀ := by
    have h := (abs_lt.mp hg4).2
    linarith
  have hsimple : TopSimple ((m.X N ω)ᵀ * m.X N ω)
      (isHermitian_transpose_mul_self (m.X N ω)) := not_not.mp g7
  -- names
  set X := m.X N ω with hXdef
  set q := m.qvec N ω with hqdef
  set vv := WithLp.ofLp (m.v N) with hvdef
  set rr := GenRMT.Deloc.rvecOf vv q with hrdef
  have hrr' : rr = vv - ((vv ⬝ᵥ q) / (q ⬝ᵥ q)) • q := rfl
  have h5 : ‖R4C.qformC (m.W0 N ω) z rr - MP.mC c z / ((m.θ : ℂ) ^ 2 + 1)‖ < δ₁ := by
    have h : ¬ (δ₁ ≤ |‖R4C.qformC (m.W0 N ω) z (vv - ((vv ⬝ᵥ q) / (q ⬝ᵥ q)) • q)
        - MP.mC c z / ((m.θ : ℂ) ^ 2 + 1)‖ - 0|) := g5
    rw [not_le, sub_zero, abs_of_nonneg (norm_nonneg _)] at h
    rw [hrr']
    exact h
  have h6 : ‖R4C.cformC (m.W0 N ω) z q rr‖ < δ₂ := by
    have h : ¬ (δ₂ ≤ |‖R4C.cformC (m.W0 N ω) z q
        (vv - ((vv ⬝ᵥ q) / (q ⬝ᵥ q)) • q)‖ - 0|) := g6
    rw [not_le, sub_zero, abs_of_nonneg (norm_nonneg _)] at h
    rw [hrr']
    exact h
  -- the decomposition, folded before the arithmetic
  have hdec := GenRMT.Deloc.overlap_le_decomp X vv q
  have hvvLp : (WithLp.toLp 2 vv : EuclideanSpace ℝ (Fin (d N))) = m.v N := rfl
  rw [hvvLp, ← hrdef] at hdec
  have hω' : ε ≤ |overlap X (m.v N) - 0| := hω
  rw [sub_zero, abs_of_nonneg (overlap_nonneg _ _)] at hω'
  -- the `q` component (the Gaussian proof, verbatim)
  have hQ2pos : 0 < R4.qform2 (m.W0 N ω) z₀ q := by
    have hhalf : (0 : ℝ) < (m.θ ^ 2 + 1) * MP.mDeriv c z₀ / 2 := by positivity
    linarith only [hhalf, hQ2]
  have hq := GenRMT.Deloc.overlap_q_le X (m.isHermitian_W0 N ω) (m.gram_eq N ω) (m.hd N)
    hsimple hlam hQ2pos
  have hovqnn : 0 ≤ overlap X (WithLp.toLp 2 q) := overlap_nonneg _ _
  have hovq : overlap X (WithLp.toLp 2 q) * R4.qform2 (m.W0 N ω) z₀ q ≤ 1 := by
    rw [le_div_iff₀ hQ2pos] at hq
    linarith only [hq]
  have ha2nn : 0 ≤ GenRMT.Deloc.aOf vv q ^ 2 := sq_nonneg _
  have hvv1 : vv ⬝ᵥ vv = 1 := m.dotProduct_v_self N
  have ha2 : GenRMT.Deloc.aOf vv q ^ 2 * (q ⬝ᵥ q) ≤ 1 := by
    have h := GenRMT.Deloc.aOf_sq_mul_le vv q
    rw [hvv1] at h
    exact h
  -- the `r` component (the window bound)
  have hrrnn : rr ⬝ᵥ rr ≤ 1 := by
    have h := GenRMT.Deloc.dotProduct_rvecOf_self_le vv q
    rw [hvv1] at h
    exact h
  have hovrnn : 0 ≤ overlap X (WithLp.toLp 2 rr) := overlap_nonneg _ _
  have hwin2 : overlap X (WithLp.toLp 2 rr)
      ≤ 2 * η * (R4C.qformC (Xᵀ * X) z rr).im := by
    refine le_trans (GenRMT.Deloc.overlap_le_specProj_Icc X hwin (WithLp.toLp 2 rr)) ?_
    exact GenRMT.Deloc.normSq_specProj_Icc_le (Xᵀ * X)
      (isHermitian_transpose_mul_self X) hη0 rr
  have him : (R4C.qformC (Xᵀ * X) z rr).im ≤ ‖R4C.qformC (Xᵀ * X) z rr‖ :=
    le_trans (le_abs_self _) (Complex.abs_im_le_norm _)
  have hid : R4C.cformC (Xᵀ * X) z rr rr
      = R4C.cformC (m.W0 N ω) z rr rr
        - R4C.cformC (Xᵀ * X) z rr q * R4C.cformC (m.W0 N ω) z q rr := by
    have h := GenRMT.Deloc.cformC_add_vecMulVec_eq (m.isHermitian_W0 N ω) hzimpos q rr rr
    rw [← m.gram_eq N ω] at h
    exact h
  have hnormq : ‖R4C.qformC (Xᵀ * X) z rr‖
      ≤ ‖R4C.qformC (m.W0 N ω) z rr‖
        + ‖R4C.cformC (Xᵀ * X) z rr q‖ * ‖R4C.cformC (m.W0 N ω) z q rr‖ := by
    have hq1 : R4C.qformC (Xᵀ * X) z rr = R4C.cformC (Xᵀ * X) z rr rr := rfl
    have hq2 : R4C.qformC (m.W0 N ω) z rr = R4C.cformC (m.W0 N ω) z rr rr := rfl
    rw [hq1, hq2, hid]
    refine le_trans (norm_sub_le _ _) ?_
    rw [norm_mul]
  have hcXq : ‖R4C.cformC (Xᵀ * X) z rr q‖ ≤ S / η := by
    have h := R4C.norm_cformC_le (isHermitian_transpose_mul_self X) hzimpos rr q
    rw [hzim] at h
    have h1 : Real.sqrt (rr ⬝ᵥ rr) ≤ 1 := by
      rw [show (1 : ℝ) = Real.sqrt 1 by simp]
      exact Real.sqrt_le_sqrt hrrnn
    have h2 : Real.sqrt (q ⬝ᵥ q) ≤ S := by
      rw [hSdef]
      exact Real.sqrt_le_sqrt hqqhi
    have h3 : Real.sqrt (rr ⬝ᵥ rr) * Real.sqrt (q ⬝ᵥ q) ≤ S :=
      le_trans (mul_le_mul h1 h2 (Real.sqrt_nonneg _) zero_le_one) (le_of_eq (one_mul S))
    refine le_trans h ?_
    exact div_le_div_of_nonneg_right h3 hη0.le
  have hqW0 : ‖R4C.qformC (m.W0 N ω) z rr‖ ≤ Cb + δ₁ := by
    have h1 : ‖R4C.qformC (m.W0 N ω) z rr‖ - ‖MP.mC c z / ((m.θ : ℂ) ^ 2 + 1)‖
        ≤ ‖R4C.qformC (m.W0 N ω) z rr - MP.mC c z / ((m.θ : ℂ) ^ 2 + 1)‖ :=
      norm_sub_norm_le _ _
    have h2 : ‖MP.mC c z / ((m.θ : ℂ) ^ 2 + 1)‖ ≤ Cb := by
      rw [norm_div]
      have hne : ‖((m.θ : ℂ) ^ 2 + 1)‖ = m.θ ^ 2 + 1 := by
        have hcast : ((m.θ : ℂ) ^ 2 + 1) = ((m.θ ^ 2 + 1 : ℝ) : ℂ) := by push_cast; ring
        rw [hcast, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hK.le]
      rw [hne, div_le_iff₀ hK]
      nlinarith only [hmCb, hK, hCb0, norm_nonneg (MP.mC c z)]
    linarith only [h1, h2, h5]
  have hstep : overlap X (WithLp.toLp 2 rr) ≤ 2 * η * ((Cb + δ₁) + (S / η) * δ₂) := by
    refine le_trans hwin2 ?_
    have hmul : (R4C.qformC (Xᵀ * X) z rr).im ≤ (Cb + δ₁) + (S / η) * δ₂ := by
      refine le_trans him (le_trans hnormq ?_)
      have hprod : ‖R4C.cformC (Xᵀ * X) z rr q‖ * ‖R4C.cformC (m.W0 N ω) z q rr‖
          ≤ (S / η) * δ₂ :=
        mul_le_mul hcXq h6.le (norm_nonneg _) (by positivity)
      linarith only [hprod, hqW0]
    have h2η : (0 : ℝ) < 2 * η := by linarith only [hη0]
    exact mul_le_mul_of_nonneg_left hmul h2η.le
  have hcancel : 2 * η * ((Cb + δ₁) + (S / η) * δ₂)
      = 2 * (η * Cb) + 2 * (η * δ₁) + 2 * (S * δ₂) := by
    field_simp
  rw [hcancel] at hstep
  -- the arithmetic. The five spectral scalars become opaque first: `linarith` and `nlinarith`
  -- compare atoms up to definitional unfolding, and `overlap` is expensive to unfold.
  set A2 : ℝ := GenRMT.Deloc.aOf vv q ^ 2 with hA2def
  set OV : ℝ := overlap X (WithLp.toLp 2 q) with hOVdef
  set OR : ℝ := overlap X (WithLp.toLp 2 rr) with hORdef
  set QF : ℝ := R4.qform2 (m.W0 N ω) z₀ q with hQFdef
  set QQ : ℝ := q ⬝ᵥ q with hQQdef
  set OL : ℝ := overlap X (m.v N) with hOLdef
  set MD : ℝ := MP.mDeriv c z₀ with hMDdef
  clear_value A2 OV OR QF QQ OL MD
  have hs1 : A2 * (m.θ ^ 2 + 1) ≤ 2 := by nlinarith only [ha2, hqqlo, ha2nn]
  have hs2 : OV * ((m.θ ^ 2 + 1) * MD) ≤ 2 := by nlinarith only [hovq, hQ2, hovqnn]
  have hs3 : (A2 * OV) * ((m.θ ^ 2 + 1) ^ 2 * MD) ≤ 4 := by
    have he : (A2 * OV) * ((m.θ ^ 2 + 1) ^ 2 * MD)
        = (A2 * (m.θ ^ 2 + 1)) * (OV * ((m.θ ^ 2 + 1) * MD)) := by ring
    rw [he]
    nlinarith only [hs1, hs2, mul_nonneg ha2nn hK.le, mul_nonneg hovqnn (mul_pos hK hmD).le]
  have hs4 : 2 * (A2 * OV) < ε / 2 := by
    have hpos : (0 : ℝ) < (m.θ ^ 2 + 1) ^ 2 * MD := by positivity
    nlinarith only [hs3, hbound, hpos, mul_nonneg ha2nn hovqnn]
  have e1 : 2 * (η * Cb) < ε / 12 := by linarith only [hηC, hη0]
  have e2 : 2 * (η * δ₁) ≤ ε / 24 := by
    rw [hδ₁def]
    nlinarith only [hη1, hη0, hε]
  have e3 : 2 * (S * δ₂) < ε / 12 := by linarith only [hδ₂e, hδ₂]
  linarith only [hdec, hω', hs4, hstep, e1, e2, e3, hε]

end SpikedModel

end StackedSVD
