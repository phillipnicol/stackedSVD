/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Deloc
import StackedSVD.RMT.General.DelocLimits
import StackedSVD.RMT.General.Defs
import StackedSVD.RMT.General.FormsBridge
import StackedSVD.RMT.General.FormsLimits
import StackedSVD.RMT.General.IsoMixed

/-!
# Item Sym at a general law: the core

The uniform delocalization bound over the unit sphere of `v ⊥`, given the limit `x` of the top
eigenvalue and a bound on the Marchenko-Pastur transform near `x`.

The route is the spectral window of `RMT/General/Deloc.lean`, as in `RMT/General/DelocAlign.lean`
(item R6'), with one difference: the direction `w` moves, so every `w`-dependent quantity needs
an explicit measure bound whose rate is free of `w` (the isotropic laws of
`RMT/General/Iso.lean` and `RMT/General/IsoMixed.lean`) rather than a `TendstoInProb` statement.

1. Window. On `|λ_max - x| ≤ η` the overlap is at most `2 η Im (wᵀ G_A(x + i η) w)`.
2. Sherman-Morrison. `A = W₀ + q qᵀ`, so `wᵀ G_A w = wᵀ G₀ w - (wᵀ G_A q)(qᵀ G₀ w)`; the first
   factor of the product is at most `√(qᵀq) / η` by the crude resolvent bound, and the second
   is small.
3. The bridge. `W₀ = gram Z - g gᵀ`, so `wᵀ G₀ w`, `wᵀ G₀ v` and `wᵀ G₀ g` are rational in
   `wᵀ G_W w`, `wᵀ G_W v`, `wᵀ G_W g` and the companion scalar `z uᵀ G̃ u`.
-/

open Filter Topology MeasureTheory ProbabilityTheory
open scoped InnerProductSpace Matrix ENNReal

namespace StackedSVD

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-! ## Two algebra helpers -/

/-- Left linearity of `cformC` at `t • a + b`. Public twin of the private helpers of
`RMT/General/Iso.lean`. -/
private theorem cformC_smul_add_left {D : ℕ} (W : Matrix (Fin D) (Fin D) ℝ) (z : ℂ) (t : ℝ)
    (a b y : Fin D → ℝ) :
    R4C.cformC W z (t • a + b) y = (t : ℂ) * R4C.cformC W z a y + R4C.cformC W z b y := by
  simp only [R4C.cformC, R4C.cvec, dotProduct, Pi.add_apply, Pi.smul_apply, smul_eq_mul,
    Complex.ofReal_add, Complex.ofReal_mul]
  rw [Finset.mul_sum, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun i _ => by ring

/-- A quotient by a quantity bounded below, with a numerator bounded above. Kept as its own
lemma so that the arithmetic in the assembly runs on small atoms. -/
private theorem div_le_mul_inv_of_le {a b t Tlo : ℝ} (hab : a ≤ b) (hb : 0 ≤ b)
    (hTlo : 0 < Tlo) (hT : Tlo ≤ t) : a / t ≤ b * Tlo⁻¹ := by
  have htpos : 0 < t := lt_of_lt_of_le hTlo hT
  have hinv : (0 : ℝ) < Tlo⁻¹ := inv_pos.mpr hTlo
  have h1 : (1 : ℝ) ≤ Tlo⁻¹ * t := by
    have h := mul_le_mul_of_nonneg_left hT hinv.le
    rwa [inv_mul_cancel₀ hTlo.ne'] at h
  rw [div_le_iff₀ htpos]
  nlinarith [mul_le_mul_of_nonneg_left h1 hb]

/-! F37 (2026-09-09): `measurable_mixDir'` used to be repeated here (`private`);
`GenRMT.measurable_mix_dir` (`IsoMixed.lean:858`, already public, already imported by this
file) is the canonical copy, used directly. -/

/-! ## The three uniform measure bounds -/

/-- (U1) The quadratic form of the noise Gram at a unit direction `w` is within `δ` of the
Stieltjes transform outside a set of measure at most `isoRate / δ²`, and the rate is free of
`w`. -/
private theorem meas_qformC_w_le (m : SpikedModel μ n d) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hG : m.GeneralNoise ν) {z : ℂ} (hz : 0 < z.im) (N : ℕ) (w : Fin (d N) → ℝ)
    (hw : w ⬝ᵥ w = 1) {δ : ℝ} (hδ : 0 < δ) :
    μ N {ω | δ ≤ ‖R4C.qformC (GenRMT.gram (m.Z N ω)) z w
        - R4C.stieltjesC (GenRMT.gram (m.Z N ω)) z‖}
      ≤ ENNReal.ofReal (GenRMT.isoRate ν z (n N) (d N) / δ ^ 2) := by
  have hmf : Measurable fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
      ‖R4C.qformC (GenRMT.gram Y) z w - R4C.stieltjesC (GenRMT.gram Y) z‖ :=
    ((GenRMT.measurable_qformC (W := fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ => GenRMT.gram Y)
      GenRMT.measurable_gram_self z (y := fun _ => w) (fun _ => measurable_const)).sub
      (GenRMT.measurable_stieltjesC GenRMT.measurable_gram_self z)).norm
  have hmset : MeasurableSet {Y : Matrix (Fin (n N)) (Fin (d N)) ℝ |
      δ ≤ ‖R4C.qformC (GenRMT.gram Y) z w - R4C.stieltjesC (GenRMT.gram Y) z‖} :=
    measurableSet_le measurable_const hmf
  have hlaw := (hG N).measure_eq (p := fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
    δ ≤ ‖R4C.qformC (GenRMT.gram Y) z w - R4C.stieltjesC (GenRMT.gram Y) z‖) hmset
  rw [show {ω : Ω N | δ ≤ ‖R4C.qformC (GenRMT.gram (m.Z N ω)) z w
      - R4C.stieltjesC (GenRMT.gram (m.Z N ω)) z‖}
      = {ω : Ω N | δ ≤ ‖R4C.qformC (GenRMT.gram (m.Z N ω)) z w
        - R4C.stieltjesC (GenRMT.gram (m.Z N ω)) z‖} from rfl, hlaw]
  have hb := GenRMT.measure_qformC_sub_ge_le hν hz (m.hn N) (m.hd N) w (le_of_eq hw) hδ
  rw [hw] at hb
  simpa using hb

/-- (U2) The bilinear form of the noise Gram at `w ⊥ v` is small outside a set of measure at
most `isoRate / δ²`. -/
private theorem meas_cformC_wv_le (m : SpikedModel μ n d) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hG : m.GeneralNoise ν) {z : ℂ} (hz : 0 < z.im) (N : ℕ) (w : Fin (d N) → ℝ)
    (hw : w ⬝ᵥ w = 1) (hwv : w ⬝ᵥ WithLp.ofLp (m.v N) = 0) {δ : ℝ} (hδ : 0 < δ) :
    μ N {ω | δ ≤ ‖R4C.cformC (GenRMT.gram (m.Z N ω)) z w (WithLp.ofLp (m.v N))‖}
      ≤ ENNReal.ofReal (GenRMT.isoRate ν z (n N) (d N) / δ ^ 2) := by
  have hmf : Measurable fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
      ‖R4C.cformC (GenRMT.gram Y) z w (WithLp.ofLp (m.v N))‖ :=
    (GenRMT.measurable_cformC (W := fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ => GenRMT.gram Y)
      GenRMT.measurable_gram_self z (x := fun _ => w)
      (y := fun _ => WithLp.ofLp (m.v N)) (fun _ => measurable_const)
      (fun _ => measurable_const)).norm
  have hmset : MeasurableSet {Y : Matrix (Fin (n N)) (Fin (d N)) ℝ |
      δ ≤ ‖R4C.cformC (GenRMT.gram Y) z w (WithLp.ofLp (m.v N))‖} :=
    measurableSet_le measurable_const hmf
  have hlaw := (hG N).measure_eq (p := fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
    δ ≤ ‖R4C.cformC (GenRMT.gram Y) z w (WithLp.ofLp (m.v N))‖) hmset
  rw [show {ω : Ω N | δ ≤ ‖R4C.cformC (GenRMT.gram (m.Z N ω)) z w (WithLp.ofLp (m.v N))‖}
      = {ω : Ω N | δ ≤ ‖R4C.cformC (GenRMT.gram (m.Z N ω)) z w (WithLp.ofLp (m.v N))‖}
      from rfl, hlaw]
  have hvv : WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) ≤ 1 := le_of_eq (m.dotProduct_v_self N)
  have hb := GenRMT.measure_cformC_sub_ge_le hν hz (m.hn N) (m.hd N) w
    (WithLp.ofLp (m.v N)) (le_of_eq hw) hvv hδ
  rw [hwv] at hb
  simpa using hb

/-- (U3) The mixed form `wᵀ G g` is small outside a set of measure at most
`isoRateMixed / δ²`. -/
private theorem meas_cformC_wg_le (m : SpikedModel μ n d) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hG : m.GeneralNoise ν) {z : ℂ} (hz : 0 < z.im) (N : ℕ) (w : Fin (d N) → ℝ)
    (hw : w ⬝ᵥ w = 1) {δ : ℝ} (hδ : 0 < δ) :
    μ N {ω | δ ≤ ‖R4C.cformC (GenRMT.gram (m.Z N ω)) z w (m.gvec N ω)‖}
      ≤ ENNReal.ofReal (GenRMT.isoRateMixed ν z (n N) (d N) / δ ^ 2) := by
  have hmf : Measurable fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
      ‖R4C.cformC (GenRMT.gram Y) z w ((Real.sqrt (d N))⁻¹ • (Yᵀ *ᵥ WithLp.ofLp (m.u N)))‖ :=
    (GenRMT.measurable_cformC (W := fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ => GenRMT.gram Y)
      GenRMT.measurable_gram_self z (x := fun _ => w)
      (y := fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
        (Real.sqrt (d N))⁻¹ • (Yᵀ *ᵥ WithLp.ofLp (m.u N)))
      (fun _ => measurable_const)
      (fun i => GenRMT.measurable_mix_dir (WithLp.ofLp (m.u N)) i)).norm
  have hmset : MeasurableSet {Y : Matrix (Fin (n N)) (Fin (d N)) ℝ |
      δ ≤ ‖R4C.cformC (GenRMT.gram Y) z w
        ((Real.sqrt (d N))⁻¹ • (Yᵀ *ᵥ WithLp.ofLp (m.u N)))‖} :=
    measurableSet_le measurable_const hmf
  have hlaw := (hG N).measure_eq (p := fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
    δ ≤ ‖R4C.cformC (GenRMT.gram Y) z w
      ((Real.sqrt (d N))⁻¹ • (Yᵀ *ᵥ WithLp.ofLp (m.u N)))‖) hmset
  have hseteq : {ω : Ω N | δ ≤ ‖R4C.cformC (GenRMT.gram (m.Z N ω)) z w (m.gvec N ω)‖}
      = {ω : Ω N | δ ≤ ‖R4C.cformC (GenRMT.gram (m.Z N ω)) z w
          ((Real.sqrt (d N))⁻¹ • ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)))‖} := by
    ext ω
    simp only [Set.mem_ofPred_eq, m.gvec_eq_smul N ω]
  rw [hseteq, hlaw]
  have huu : WithLp.ofLp (m.u N) ⬝ᵥ WithLp.ofLp (m.u N) ≤ 1 := le_of_eq (m.dotProduct_u_self N)
  exact GenRMT.measure_cformC_mixed_ge_le hν hz (m.hn N) (m.hd N) w (le_of_eq hw) _ huu hδ

/-! ## The theorem -/

/-- The arithmetic of the assembly, on plain reals: the window height `η`, the two tolerances
`δ` and `σ`, and the four constants `L`, `Sq`, `θ`, `B` give a bound below `ε`. Kept apart so
that `nlinarith` never sees a spectral atom. -/
private theorem window_lt {ε η mCn δ σ L Sq θ B t : ℝ}
    (_hε : 0 < ε) (_hη0 : 0 < η) (hη1 : η ≤ 1) (hmC : η * mCn < ε / 8)
    (hδ0 : 0 < δ) (hδ1 : δ ≤ 1) (hσ : σ = ε / 24) (hδs : δ ≤ ε / 100)
    (hθ : 0 ≤ θ) (_hSq0 : 0 ≤ Sq) (hL0 : 0 < L) (hB1 : 1 ≤ B)
    (hBS : Sq ≤ B) (hBθ : θ ≤ B) (hBL : L ≤ B) (hBδ : B * δ ≤ ε / 100)
    (hBBδ : B * B * δ ≤ ε / 100)
    (ht : t ≤ 2 * η * mCn + 2 * η * (δ + σ + δ * δ * L) + 2 * (Sq * (θ * (2 * δ) + δ * L))) :
    t < ε := by
  have hB0 : (0 : ℝ) < B := by linarith
  have hd1 : δ * δ * L ≤ δ * δ * B := by nlinarith [mul_nonneg hδ0.le hδ0.le]
  have hd2 : δ * δ * B ≤ B * δ := by
    nlinarith [mul_nonneg (mul_nonneg hB0.le hδ0.le) (by linarith : (0 : ℝ) ≤ 1 - δ)]
  have hdd : δ * δ * L ≤ ε / 100 := by linarith
  have hp1 : 2 * η * mCn < ε / 4 := by linarith
  have hσ0 : (0 : ℝ) ≤ σ := by rw [hσ]; linarith
  have hnn : (0 : ℝ) ≤ δ + σ + δ * δ * L := by
    have h : (0 : ℝ) ≤ δ * δ * L := by positivity
    linarith
  have hle : δ + σ + δ * δ * L ≤ ε / 8 := by rw [hσ]; linarith
  have hp2 : 2 * η * (δ + σ + δ * δ * L) ≤ ε / 4 := by
    nlinarith [mul_le_mul_of_nonneg_right hη1 hnn]
  have hbr : θ * (2 * δ) + δ * L ≤ 3 * (B * δ) := by
    nlinarith [mul_le_mul_of_nonneg_right hBθ (by linarith : (0 : ℝ) ≤ 2 * δ),
      mul_le_mul_of_nonneg_left hBL hδ0.le]
  have hbrnn : (0 : ℝ) ≤ θ * (2 * δ) + δ * L := by positivity
  have hstep3 : Sq * (θ * (2 * δ) + δ * L) ≤ B * (3 * (B * δ)) :=
    mul_le_mul hBS hbr hbrnn (by linarith)
  have hp3 : 2 * (Sq * (θ * (2 * δ) + δ * L)) ≤ ε / 4 := by linarith
  linarith

set_option maxHeartbeats 1000000 in
-- The assembly carries five bad families, three moving families and eight spectral scalars in
-- one block, so the default budget is not enough.
/-- **Item Sym at a general law, the core.** Uniform over the unit sphere of `v ⊥`, given the
limit `x` of the top eigenvalue and a bound on the MP transform near `x`. -/
theorem delocUniform_of_lamMax [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (H : m.ResolventFormsC c) (hreg : m.Regime c) {ν : Measure ℝ}
    (hν : NoiseLaw ν) (hG : m.GeneralNoise ν) {x : ℝ}
    (hlam : TendstoInProb μ (fun N ω => gramLamMax (m.X N ω)) x)
    (hm : ∀ δ > 0, ∃ η : ℝ, 0 < η ∧ η ≤ 1 ∧
      η * ‖MP.mC c ((x : ℂ) + (η : ℂ) * Complex.I)‖ < δ) :
    ∀ ε > 0, Tendsto (fun N => ⨆ w ∈ m.orthUnit N, μ N {ω | ε ≤ overlap (m.X N ω) w})
      atTop (𝓝 0) := by
  intro ε hε
  obtain ⟨η, hη0, hη1, hηm⟩ := hm (ε / 8) (by positivity)
  set z : ℂ := (x : ℂ) + (η : ℂ) * Complex.I with hzdef
  have hzim : z.im = η := by simp [hzdef]
  have hz : 0 < z.im := by rw [hzim]; exact hη0
  have hzne : z ≠ 0 := MP.ne_zero_of_im_pos hz
  have hzn : 0 < ‖z‖ := norm_pos_iff.mpr hzne
  -- the nonzero limit of the companion denominator
  have hT0ne : z * MP.mTildeC c z ≠ 0 := MP.z_mul_mTildeC_ne_zero hc hz
  set Tlo : ℝ := ‖z * MP.mTildeC c z‖ / 2 with hTlodef
  have hTlo : 0 < Tlo := by
    rw [hTlodef]
    exact half_pos (norm_pos_iff.mpr hT0ne)
  set L : ℝ := Tlo⁻¹ with hLdef
  have hL0 : 0 < L := by rw [hLdef]; exact inv_pos.mpr hTlo
  have hLT : L * Tlo = 1 := by rw [hLdef]; field_simp
  set Sq : ℝ := Real.sqrt (m.θ ^ 2 + 2) with hSqdef
  have hSq0 : 0 ≤ Sq := Real.sqrt_nonneg _
  set B : ℝ := 1 + m.θ + Sq + L with hBdef
  have hB1 : 1 ≤ B := by rw [hBdef]; linarith [m.hθ, hSq0, hL0.le]
  have hB0 : 0 < B := by linarith
  have hBθ : m.θ ≤ B := by rw [hBdef]; linarith [hSq0, hL0.le]
  have hBS : Sq ≤ B := by rw [hBdef]; linarith [m.hθ, hL0.le]
  have hBL : L ≤ B := by rw [hBdef]; linarith [m.hθ, hSq0]
  set δ : ℝ := min 1 (ε / (100 * B * B)) with hδdef
  have hδ0 : 0 < δ := lt_min one_pos (by positivity)
  have hδ1 : δ ≤ 1 := min_le_left _ _
  have hδB : 100 * B * B * δ ≤ ε := by
    have h : δ ≤ ε / (100 * B * B) := min_le_right _ _
    rw [le_div_iff₀ (by positivity)] at h
    linarith
  have hBBδ : B * B * δ ≤ ε / 100 := by linarith
  have hBδ : B * δ ≤ ε / 100 := by nlinarith [hB1, hδ0.le, hB0]
  have hδs : δ ≤ ε / 100 := by nlinarith [hB1, hδ0.le]
  set σ : ℝ := ε / 24 with hσdef
  have hσ0 : 0 < σ := by rw [hσdef]; positivity
  -- the five `w`-independent bad families
  set Bad : ∀ N, Set (Ω N) := fun N =>
    {ω | η ≤ |gramLamMax (m.X N ω) - x|}
      ∪ ({ω | σ ≤ ‖R4C.stieltjesC (GenRMT.gram (m.Z N ω)) z - MP.mC c z‖}
      ∪ ({ω | Tlo / ‖z‖ ≤ ‖R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))
            - MP.mTildeC c z‖}
      ∪ ({ω | (1 : ℝ) ≤ ‖R4C.cformC (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω)‖}
      ∪ {ω | (1 : ℝ) ≤ |m.qvec N ω ⬝ᵥ m.qvec N ω - (m.θ ^ 2 + 1)|}))) with hBaddef
  have hb1 : Tendsto (fun N => μ N {ω : Ω N | η ≤ |gramLamMax (m.X N ω) - x|}) atTop (𝓝 0) :=
    hlam η hη0
  have hb2 : Tendsto (fun N => μ N {ω : Ω N |
      σ ≤ ‖R4C.stieltjesC (GenRMT.gram (m.Z N ω)) z - MP.mC c z‖}) atTop (𝓝 0) := by
    have h := GenRMT.tendstoInProb_stieltjesC_general (μ := μ) hν hc hz hreg.2.1 m.hn hreg.2.2
      m.Z hG σ hσ0
    refine h.congr fun N => ?_
    congr 1
    ext ω
    simp only [Set.mem_ofPred_eq, sub_zero, abs_of_nonneg (norm_nonneg _)]
  have hb3 : Tendsto (fun N => μ N {ω : Ω N | Tlo / ‖z‖
      ≤ ‖R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)) - MP.mTildeC c z‖})
      atTop (𝓝 0) := by
    have h := m.qformC_gramC_u_tendsto hc hreg hν hG z hz (Tlo / ‖z‖) (by positivity)
    refine h.congr fun N => ?_
    congr 1
    ext ω
    simp only [Set.mem_ofPred_eq, sub_zero, abs_of_nonneg (norm_nonneg _)]
  have hb4 : Tendsto (fun N => μ N {ω : Ω N |
      (1 : ℝ) ≤ ‖R4C.cformC (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω)‖})
      atTop (𝓝 0) := by
    have h := H.vgC z hz 1 one_pos
    refine h.congr fun N => ?_
    congr 1
    ext ω
    simp only [Set.mem_ofPred_eq, sub_zero, abs_of_nonneg (norm_nonneg _)]
  have hb5 : Tendsto (fun N => μ N {ω : Ω N |
      (1 : ℝ) ≤ |m.qvec N ω ⬝ᵥ m.qvec N ω - (m.θ ^ 2 + 1)|}) atTop (𝓝 0) :=
    m.tendstoInProb_dotProduct_qvec_general hν hG hreg.2.1 1 one_pos
  have hBad : Tendsto (fun N => μ N (Bad N)) atTop (𝓝 0) :=
    tendsto_measure_zero_union hb1 (tendsto_measure_zero_union hb2
      (tendsto_measure_zero_union hb3 (tendsto_measure_zero_union hb4 hb5)))
  -- the rate of the three `w`-dependent families
  set rate : ℕ → ℝ≥0∞ := fun N =>
    ENNReal.ofReal (GenRMT.isoRate ν z (n N) (d N) / δ ^ 2)
      + (ENNReal.ofReal (GenRMT.isoRate ν z (n N) (d N) / δ ^ 2)
      + ENNReal.ofReal (GenRMT.isoRateMixed ν z (n N) (d N) / δ ^ 2)) with hratedef
  have hr1 : Tendsto (fun N => ENNReal.ofReal (GenRMT.isoRate ν z (n N) (d N) / δ ^ 2)) atTop
      (𝓝 0) := by
    have h : Tendsto (fun N => GenRMT.isoRate ν z (n N) (d N) / δ ^ 2) atTop (𝓝 0) := by
      have h0 := (GenRMT.tendsto_isoRate ν hz hreg.2.1 hreg.2.2).div_const (δ ^ 2)
      simpa using h0
    have h2 := ENNReal.tendsto_ofReal h
    simpa using h2
  have hr2 : Tendsto (fun N => ENNReal.ofReal (GenRMT.isoRateMixed ν z (n N) (d N) / δ ^ 2))
      atTop (𝓝 0) := by
    have h : Tendsto (fun N => GenRMT.isoRateMixed ν z (n N) (d N) / δ ^ 2) atTop (𝓝 0) := by
      have h0 := (GenRMT.tendsto_isoRateMixed ν hz hreg.2.1 hreg.2.2).div_const (δ ^ 2)
      simpa using h0
    have h2 := ENNReal.tendsto_ofReal h
    simpa using h2
  have hrate : Tendsto rate atTop (𝓝 0) := by
    rw [hratedef]
    have h := hr1.add (hr1.add hr2)
    simpa using h
  have htot : Tendsto (fun N => μ N (Bad N) + rate N) atTop (𝓝 0) := by
    have h := hBad.add hrate
    simpa using h
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds htot
    (fun _ => zero_le) fun N => ?_
  refine iSup₂_le fun w hw => ?_
  -- the two facts about `w`
  obtain ⟨hw1, hwv⟩ := hw
  have hww : WithLp.ofLp w ⬝ᵥ WithLp.ofLp w = 1 := by
    have h : ⟪w, w⟫_ℝ = WithLp.ofLp w ⬝ᵥ WithLp.ofLp w := inner_euclidean_eq_dotProduct w w
    rw [real_inner_self_eq_norm_sq, hw1] at h
    simpa using h.symm
  have hwvd : WithLp.ofLp w ⬝ᵥ WithLp.ofLp (m.v N) = 0 := by
    have h : ⟪w, m.v N⟫_ℝ = WithLp.ofLp w ⬝ᵥ WithLp.ofLp (m.v N) :=
      inner_euclidean_eq_dotProduct w (m.v N)
    rw [hwv] at h
    exact h.symm
  -- the inclusion
  have hsub : {ω : Ω N | ε ≤ overlap (m.X N ω) w}
      ⊆ Bad N
        ∪ ({ω | δ ≤ ‖R4C.qformC (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp w)
              - R4C.stieltjesC (GenRMT.gram (m.Z N ω)) z‖}
        ∪ ({ω | δ ≤ ‖R4C.cformC (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp w)
              (WithLp.ofLp (m.v N))‖}
        ∪ {ω | δ ≤ ‖R4C.cformC (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp w) (m.gvec N ω)‖})) := by
    intro ω hω
    by_contra hcon
    simp only [Set.mem_union, not_or] at hcon
    obtain ⟨gBad, c1, c2, c3⟩ := hcon
    rw [hBaddef] at gBad
    simp only [Set.mem_union, not_or, Set.mem_ofPred_eq] at gBad
    obtain ⟨g1, g2, g3, g4, g5⟩ := gBad
    rw [not_le] at g1 g2 g3 g4 g5
    rw [Set.mem_ofPred_eq, not_le] at c1 c2 c3
    -- names
    set W : Matrix (Fin (d N)) (Fin (d N)) ℝ := GenRMT.gram (m.Z N ω) with hWdef
    have hWh : W.IsHermitian := GenRMT.gram_isHermitian (m.Z N ω)
    have hW0h : (m.W0 N ω).IsHermitian := m.isHermitian_W0 N ω
    set g : Fin (d N) → ℝ := m.gvec N ω with hgdef
    set q : Fin (d N) → ℝ := m.qvec N ω with hqdef
    set vv : Fin (d N) → ℝ := WithLp.ofLp (m.v N) with hvvdef
    set ww : Fin (d N) → ℝ := WithLp.ofLp w with hwwdef
    set T : ℂ := z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)) with hTdef
    -- the denominator is bounded below
    have hTsub : ‖T - z * MP.mTildeC c z‖ < Tlo := by
      have he : T - z * MP.mTildeC c z
          = z * (R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))
            - MP.mTildeC c z) := by rw [hTdef]; ring
      rw [he, norm_mul]
      calc ‖z‖ * ‖R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))
            - MP.mTildeC c z‖ < ‖z‖ * (Tlo / ‖z‖) := by
            exact mul_lt_mul_of_pos_left g3 hzn
        _ = Tlo := by field_simp
    have hTn : Tlo ≤ ‖T‖ := by
      have h1 : ‖z * MP.mTildeC c z‖ - ‖T‖ ≤ ‖z * MP.mTildeC c z - T‖ := norm_sub_norm_le _ _
      rw [norm_sub_rev] at h1
      have h2 : ‖z * MP.mTildeC c z‖ = 2 * Tlo := by rw [hTlodef]; ring
      linarith
    have hTpos : 0 < ‖T‖ := lt_of_lt_of_le hTlo hTn
    have hTne : T ≠ 0 := norm_pos_iff.mp hTpos
    -- (1) the quadratic form on `W₀` at `w`
    have hq0 : R4C.qformC (m.W0 N ω) z ww
        = R4C.qformC W z ww - R4C.cformC W z ww g * R4C.cformC W z ww g / T := by
      have h := m.qformC_W0_eq N ω hz ww
      rw [GenRMT.cformC_comm hWh hz.ne' g ww] at h
      exact h
    have hnum : ‖R4C.cformC W z ww g * R4C.cformC W z ww g / T‖ ≤ δ * δ * L := by
      rw [norm_div, norm_mul, hLdef]
      refine div_le_mul_inv_of_le ?_ (by positivity) hTlo hTn
      exact mul_le_mul c3.le c3.le (norm_nonneg _) hδ0.le
    have hqW0 : ‖R4C.qformC (m.W0 N ω) z ww - MP.mC c z‖ ≤ δ + σ + δ * δ * L := by
      have he : R4C.qformC (m.W0 N ω) z ww - MP.mC c z
          = (R4C.qformC W z ww - R4C.stieltjesC W z)
            + (R4C.stieltjesC W z - MP.mC c z)
            - R4C.cformC W z ww g * R4C.cformC W z ww g / T := by
        rw [hq0]; ring
      rw [he]
      refine le_trans (norm_sub_le _ _) ?_
      refine add_le_add (le_trans (norm_add_le _ _) ?_) hnum
      exact add_le_add c1.le g2.le
    -- (2) the mixed form on `W₀`
    have hcg : R4C.cformC (m.W0 N ω) z ww g = -(R4C.cformC W z ww g / T) :=
      m.cformC_W0_gvec_eq N ω hz ww
    have hcgb : ‖R4C.cformC (m.W0 N ω) z ww g‖ ≤ δ * L := by
      rw [hcg, norm_neg, norm_div, hLdef]
      exact div_le_mul_inv_of_le c3.le hδ0.le hTlo hTn
    -- (3) the form on `W₀` at `v`
    have hgram : m.W0 N ω + Matrix.vecMulVec g g = W := by
      rw [hgdef, hWdef, m.W0_eq_gram_sub N ω]
      abel
    have hcv : R4C.cformC (m.W0 N ω) z ww vv
        = R4C.cformC W z ww vv + R4C.cformC W z ww g * R4C.cformC (m.W0 N ω) z g vv := by
      have h := GenRMT.Deloc.cformC_add_vecMulVec_eq hW0h hz g ww vv
      rw [hgram] at h
      linear_combination -h
    have hcvb : ‖R4C.cformC (m.W0 N ω) z ww vv‖ ≤ 2 * δ := by
      rw [hcv]
      refine le_trans (norm_add_le _ _) ?_
      have h2 : ‖R4C.cformC W z ww g * R4C.cformC (m.W0 N ω) z g vv‖ ≤ δ * 1 := by
        rw [norm_mul]
        have hgv : ‖R4C.cformC (m.W0 N ω) z g vv‖ < 1 := by
          rw [GenRMT.cformC_comm hW0h hz.ne' g vv]
          exact g4
        exact mul_le_mul c3.le hgv.le (norm_nonneg _) hδ0.le
      linarith [c2.le, h2]
    -- (4) the form on `W₀` at `q`
    have hqw : R4C.cformC (m.W0 N ω) z q ww
        = (m.θ : ℂ) * R4C.cformC (m.W0 N ω) z ww vv + R4C.cformC (m.W0 N ω) z ww g := by
      have hqe : q = m.θ • vv + g := rfl
      rw [hqe, cformC_smul_add_left]
      rw [GenRMT.cformC_comm hW0h hz.ne' vv ww, GenRMT.cformC_comm hW0h hz.ne' g ww]
    have hqwb : ‖R4C.cformC (m.W0 N ω) z q ww‖ ≤ m.θ * (2 * δ) + δ * L := by
      rw [hqw]
      refine le_trans (norm_add_le _ _) ?_
      have h1 : ‖(m.θ : ℂ) * R4C.cformC (m.W0 N ω) z ww vv‖ ≤ m.θ * (2 * δ) := by
        rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg m.hθ]
        exact mul_le_mul_of_nonneg_left hcvb m.hθ
      linarith [hcgb]
    -- (5) Sherman-Morrison at `A`
    have hqq : q ⬝ᵥ q ≤ m.θ ^ 2 + 2 := by
      have h := (abs_lt.mp g5).2
      linarith
    have hAh : ((m.X N ω)ᵀ * m.X N ω).IsHermitian := isHermitian_transpose_mul_self (m.X N ω)
    have hsm : R4C.qformC ((m.X N ω)ᵀ * m.X N ω) z ww
        = R4C.qformC (m.W0 N ω) z ww
          - R4C.cformC ((m.X N ω)ᵀ * m.X N ω) z ww q * R4C.cformC (m.W0 N ω) z q ww := by
      have h := GenRMT.Deloc.cformC_add_vecMulVec_eq hW0h hz q ww ww
      rw [← m.gram_eq N ω] at h
      exact h
    have hAq : ‖R4C.cformC ((m.X N ω)ᵀ * m.X N ω) z ww q‖ ≤ Sq / η := by
      have h := R4C.norm_cformC_le hAh hz ww q
      rw [hzim, hww, Real.sqrt_one, one_mul] at h
      have h2 : Real.sqrt (q ⬝ᵥ q) ≤ Sq := by rw [hSqdef]; exact Real.sqrt_le_sqrt hqq
      exact le_trans h (div_le_div_of_nonneg_right h2 hη0.le)
    have hAn : ‖R4C.qformC ((m.X N ω)ᵀ * m.X N ω) z ww‖
        ≤ (δ + σ + δ * δ * L + ‖MP.mC c z‖) + (Sq / η) * (m.θ * (2 * δ) + δ * L) := by
      rw [hsm]
      refine le_trans (norm_sub_le _ _) (add_le_add ?_ ?_)
      · have h1 : ‖R4C.qformC (m.W0 N ω) z ww‖ - ‖MP.mC c z‖
            ≤ ‖R4C.qformC (m.W0 N ω) z ww - MP.mC c z‖ := norm_sub_norm_le _ _
        linarith [hqW0]
      · rw [norm_mul]
        exact mul_le_mul hAq hqwb (norm_nonneg _) (by positivity)
    -- (6) the window
    have hwin : |gramLamMax (m.X N ω) - x| ≤ η := g1.le
    have hovl : overlap (m.X N ω) w
        ≤ 2 * η * (R4C.qformC ((m.X N ω)ᵀ * m.X N ω) z ww).im := by
      refine le_trans (GenRMT.Deloc.overlap_le_specProj_Icc (m.X N ω) hwin w) ?_
      have h := GenRMT.Deloc.normSq_specProj_Icc_le ((m.X N ω)ᵀ * m.X N ω) hAh (x := x) hη0 ww
      have hlp : (WithLp.toLp 2 ww : EuclideanSpace ℝ (Fin (d N))) = w := rfl
      rw [hlp] at h
      exact h
    have him : (R4C.qformC ((m.X N ω)ᵀ * m.X N ω) z ww).im
        ≤ ‖R4C.qformC ((m.X N ω)ᵀ * m.X N ω) z ww‖ :=
      le_trans (le_abs_self _) (Complex.abs_im_le_norm _)
    -- the arithmetic
    have hεov : ε ≤ overlap (m.X N ω) w := hω
    have hstep : overlap (m.X N ω) w
        ≤ 2 * η * ((δ + σ + δ * δ * L + ‖MP.mC c z‖) + (Sq / η) * (m.θ * (2 * δ) + δ * L)) := by
      refine le_trans hovl ?_
      exact mul_le_mul_of_nonneg_left (le_trans him hAn) (by positivity)
    have hcancel : 2 * η * ((δ + σ + δ * δ * L + ‖MP.mC c z‖)
          + (Sq / η) * (m.θ * (2 * δ) + δ * L))
        = 2 * η * ‖MP.mC c z‖ + 2 * η * (δ + σ + δ * δ * L)
          + 2 * (Sq * (m.θ * (2 * δ) + δ * L)) := by
      field_simp
      ring
    rw [hcancel] at hstep
    -- the arithmetic, on plain reals
    have hfin := window_lt hε hη0 hη1 hηm hδ0 hδ1 hσdef hδs m.hθ hSq0 hL0 hB1 hBS hBθ hBL
      hBδ hBBδ hstep
    linarith [hεov, hfin]
  -- pass to the bound
  refine le_trans (measure_mono hsub) ?_
  refine le_trans (measure_union_le _ _) ?_
  refine add_le_add le_rfl ?_
  refine le_trans (measure_union_le _ _) ?_
  simp only [hratedef]
  refine add_le_add (meas_qformC_w_le m hν hG hz N (WithLp.ofLp w) hww hδ0) ?_
  refine le_trans (measure_union_le _ _) ?_
  exact add_le_add (meas_cformC_wv_le m hν hG hz N (WithLp.ofLp w) hww hwvd hδ0)
    (meas_cformC_wg_le m hν hG hz N (WithLp.ofLp w) hww hδ0)

end SpikedModel

end StackedSVD
