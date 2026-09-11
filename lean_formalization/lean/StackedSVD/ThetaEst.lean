/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Full
import StackedSVD.SVDStack.Gram
import StackedSVD.Prob.NoiseMoments

/-!
# `thm:theta_est`: estimating a signal strength below the detection threshold

`main_paper.tex:958` (statement), `main_paper.tex:2442` (proof),
`main_paper.tex:2496` (item P, `lem:noise_projection_concentration`).
Review note: `notes/archive/thm_theta_est.md`.

Table `i` is above its threshold (`c_i < θ_i⁴`). Its top singular value gives a consistent
`θ̂_i` by inverting `σ₁²(X_i) → θ_i² + 1 + c_i + c_i/θ_i²`, hence a consistent `β̂_i`. The
projection of table `j` on `v̂_i` then estimates `θ_j`, whatever the size of `θ_j`:

  `θ̂_j = √(‖X_j v̂_i‖² - c_j) / β̂_i`.

## Content

1. Scalars: `thetaHat1` (`eq:theta_est_quadratic`), `betaHat`, `thetaHatFn`
   (`eq:theta_estimation` as a function of the two observables), the inversion
   `thetaHat1_rhoSq`, and the continuity of `thetaHatFn` at the limit point.
2. `topDir X`, a measurable unit top right singular vector of `X`. It is `delocDir 0 X` of
   `SVDStack/Gram.lean`, since `perpOf 0` is the identity; the project has no measurability
   statement for `vMax`.
3. `measurePreserving_mulVec_unit`: for a unit `a`, `Z ↦ Z a` carries `gaussianMatrix r D` to
   the standard Gaussian vector law of `ℝ^r`. Two second moments follow, uniform in `a`:
   `∫ (d⁻¹ ‖Z a‖² - r/d)² = varSq · r / d²` and `∫ (x ⬝ᵥ Z a)² = 1`.
4. `SpikedModel.noise_projection_tendsto`: item P for a deterministic direction.
5. `MultiTableModel.measure_randomDir_le`: the Fubini step of item P for a direction that
   reads another table, through `pi_pair_le`, the route of `lem_delocalization`.
6. `MultiTableModel.thetaHat`, `thm_theta_est` and the Gaussian corollary.

The variance constant is `R2.varSq = ∫ (t²-1)² dN(0,1)` (its value `2` is never used), the
same opaque constant that `RMT/R2.lean` uses.

STATUS: see `notes/archive/agent_reports/theta_est.md`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal

namespace StackedSVD

/-! ### Two closure lemmas for `TendstoInProb`

`Prob/TendstoInProb.lean` has no lemma that adds a convergent deterministic sequence, and no
lemma that reads a Chebyshev bound. Both are needed here because the natural center of
`‖E a‖²` is `n N / d N`, not the limit `c`. -/

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}

/-- A vanishing random part plus a convergent deterministic part. -/
theorem TendstoInProb.add_tendsto {g : ∀ N, Ω N → ℝ} {h : ℕ → ℝ} {b : ℝ}
    (hg : TendstoInProb μ g 0) (hh : Tendsto h atTop (𝓝 b)) :
    TendstoInProb μ (fun N ω => g N ω + h N) b := by
  intro ε hε
  obtain ⟨N₀, hN₀⟩ := Metric.tendsto_atTop.mp hh (ε / 2) (by positivity)
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds
    (hg (ε / 2) (by positivity)) (Eventually.of_forall fun _ => zero_le) ?_
  filter_upwards [eventually_ge_atTop N₀] with N hN
  refine measure_mono fun ω hω => ?_
  have h1 : ε ≤ |g N ω + h N - b| := hω
  have h2 : |h N - b| < ε / 2 := by
    have h3 := hN₀ N hN
    rwa [Real.dist_eq] at h3
  change ε / 2 ≤ |g N ω - 0|
  rw [sub_zero]
  have h4 : |g N ω + h N - b| ≤ |g N ω| + |h N - b| := by
    have he : g N ω + h N - b = g N ω + (h N - b) := by ring
    rw [he]
    exact abs_add_le _ _
  linarith

namespace ThetaEst

/-- Chebyshev bounds with a vanishing constant give a limit in probability. -/
theorem tendstoInProb_of_meas_le {f : ∀ N, Ω N → ℝ} {ctr : ℕ → ℝ} {K : ℕ → ℝ}
    (hK : Tendsto K atTop (𝓝 0))
    (hb : ∀ ε : ℝ, 0 < ε → ∀ N, μ N {ω | ε ≤ |f N ω - ctr N|}
      ≤ ENNReal.ofReal (K N / ε ^ 2)) :
    TendstoInProb μ (fun N ω => f N ω - ctr N) 0 := by
  intro ε hε
  have h1 : Tendsto (fun N => ENNReal.ofReal (K N / ε ^ 2)) atTop (𝓝 0) := by
    have h2 : Tendsto (fun N => K N / ε ^ 2) atTop (𝓝 0) := by
      simpa using hK.div_const (ε ^ 2)
    simpa using ENNReal.tendsto_ofReal h2
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds h1
    (fun _ => zero_le) fun N => ?_
  have hset : {ω : Ω N | ε ≤ |f N ω - ctr N - 0|} = {ω | ε ≤ |f N ω - ctr N|} := by
    simp
  rw [hset]
  exact hb ε hε N

/-! ### The scalars of `eq:theta_estimation` -/

/-- `eq:theta_est_quadratic`: the larger root of `t² - (s - (1+c)) t + c = 0`, then its square
root. At `s = ρ²(θ, c)` and above the threshold it returns `θ` (`thetaHat1_rhoSq`). -/
noncomputable def thetaHat1 (s c : ℝ) : ℝ :=
  Real.sqrt ((s - (1 + c) + Real.sqrt ((s - (1 + c)) ^ 2 - 4 * c)) / 2)

/-- `β̂₁`, the plug-in estimate of `β₁ = √(betaSq θ₁ c₁)`. -/
noncomputable def betaHat (s c : ℝ) : ℝ := Real.sqrt (betaSq (thetaHat1 s c) c)

/-- `eq:theta_estimation`, as a function of the two observables `s = σ₁²(X_i)` and
`p = ‖X_j v̂_i‖²`. -/
noncomputable def thetaHatFn (ci cj : ℝ) (sp : ℝ × ℝ) : ℝ :=
  Real.sqrt (sp.2 - cj) / betaHat sp.1 ci

theorem continuous_thetaHat1 (c : ℝ) : Continuous fun s => thetaHat1 s c := by
  unfold thetaHat1
  refine Real.continuous_sqrt.comp (Continuous.div_const ?_ 2)
  refine (continuous_id.sub continuous_const).add (Real.continuous_sqrt.comp ?_)
  exact ((continuous_id.sub continuous_const).pow 2).sub continuous_const

/-- **The inversion.** `thetaHat1` undoes `rhoSq` above the threshold. Checked with sympy in a
session script (`check_inversion.py`, not kept). -/
theorem thetaHat1_rhoSq {θ c : ℝ} (hθ : 0 ≤ θ) (hc : 0 ≤ c) (hthr : c < θ ^ 4) :
    thetaHat1 (rhoSq θ c) c = θ := by
  have h4 : 0 < θ ^ 4 := lt_of_le_of_lt hc hthr
  have h2 : 0 < θ ^ 2 := by nlinarith [sq_nonneg θ, sq_nonneg (θ ^ 2)]
  have hθ0 : θ ≠ 0 := by
    intro h
    rw [h] at h2
    simp at h2
  have h1 : θ ^ 2 + 1 + c + c / θ ^ 2 - (1 + c) = θ ^ 2 + c / θ ^ 2 := by ring
  have hdisc : (θ ^ 2 + c / θ ^ 2) ^ 2 - 4 * c = (θ ^ 2 - c / θ ^ 2) ^ 2 := by
    field_simp
    ring
  have hpos : 0 ≤ θ ^ 2 - c / θ ^ 2 := by
    rw [sub_nonneg, div_le_iff₀ h2]
    nlinarith
  rw [rhoSq, if_pos hthr, thetaHat1, h1, hdisc, Real.sqrt_sq hpos]
  have hhalf : (θ ^ 2 + c / θ ^ 2 + (θ ^ 2 - c / θ ^ 2)) / 2 = θ ^ 2 := by ring
  rw [hhalf, Real.sqrt_sq hθ]

/-- `betaSq` is positive above the threshold. -/
theorem betaSq_pos {θ c : ℝ} (hc : 0 ≤ c) (hthr : c < θ ^ 4) : 0 < betaSq θ c := by
  have h4 : 0 < θ ^ 4 := lt_of_le_of_lt hc hthr
  have h2 : 0 < θ ^ 2 := by nlinarith [sq_nonneg θ, sq_nonneg (θ ^ 2)]
  rw [betaSq, if_pos hthr]
  exact div_pos (by linarith) (by linarith)

/-- `betaSq · c` is continuous at every point above the threshold: there the branch test is
locally constant. -/
theorem continuousAt_betaSq {θ c : ℝ} (hc : 0 ≤ c) (hthr : c < θ ^ 4) :
    ContinuousAt (fun x => betaSq x c) θ := by
  have h4 : 0 < θ ^ 4 := lt_of_le_of_lt hc hthr
  have h2 : 0 < θ ^ 2 := by nlinarith [sq_nonneg θ, sq_nonneg (θ ^ 2)]
  have hden : θ ^ 4 + θ ^ 2 ≠ 0 := by
    intro h
    linarith
  have hopen : IsOpen {x : ℝ | c < x ^ 4} :=
    isOpen_lt continuous_const (continuous_id.pow 4)
  have hU : {x : ℝ | c < x ^ 4} ∈ 𝓝 θ := hopen.mem_nhds hthr
  refine ContinuousAt.congr (f := fun x : ℝ => (x ^ 4 - c) / (x ^ 4 + x ^ 2)) ?_ ?_
  · exact ContinuousAt.div (by fun_prop) (by fun_prop) hden
  · filter_upwards [hU] with x hx
    rw [betaSq, if_pos hx]

theorem betaHat_rhoSq {θ c : ℝ} (hθ : 0 ≤ θ) (hc : 0 ≤ c) (hthr : c < θ ^ 4) :
    betaHat (rhoSq θ c) c = Real.sqrt (betaSq θ c) := by
  rw [betaHat, thetaHat1_rhoSq hθ hc hthr]

theorem continuousAt_betaHat {θ c : ℝ} (hθ : 0 ≤ θ) (hc : 0 ≤ c) (hthr : c < θ ^ 4) :
    ContinuousAt (fun s => betaHat s c) (rhoSq θ c) := by
  have h1 : ContinuousAt (fun s => thetaHat1 s c) (rhoSq θ c) :=
    (continuous_thetaHat1 c).continuousAt
  have h2 : ContinuousAt (fun x => betaSq x c) (thetaHat1 (rhoSq θ c) c) := by
    rw [thetaHat1_rhoSq hθ hc hthr]
    exact continuousAt_betaSq hc hthr
  have h3 : ContinuousAt (fun s => betaSq (thetaHat1 s c) c) (rhoSq θ c) :=
    ContinuousAt.comp (g := fun x => betaSq x c) (f := fun s => thetaHat1 s c)
      (x := rhoSq θ c) h2 h1
  have h4 := ContinuousAt.comp (g := Real.sqrt)
    (f := fun s => betaSq (thetaHat1 s c) c) (x := rhoSq θ c)
    Real.continuous_sqrt.continuousAt h3
  simpa [betaHat, Function.comp_def] using h4

/-- The value of the estimator at the limit point: `θ_j`. -/
theorem thetaHatFn_limit {θi θj ci cj : ℝ} (hθi : 0 ≤ θi) (hθj : 0 ≤ θj) (hci : 0 ≤ ci)
    (hthr : ci < θi ^ 4) :
    thetaHatFn ci cj (rhoSq θi ci, θj ^ 2 * betaSq θi ci + cj) = θj := by
  have hb : 0 < betaSq θi ci := betaSq_pos hci hthr
  have hsb : 0 < Real.sqrt (betaSq θi ci) := Real.sqrt_pos.mpr hb
  change Real.sqrt (θj ^ 2 * betaSq θi ci + cj - cj) / betaHat (rhoSq θi ci) ci = θj
  rw [betaHat_rhoSq hθi hci hthr, add_sub_cancel_right, Real.sqrt_mul (sq_nonneg θj),
    Real.sqrt_sq hθj, mul_div_assoc, div_self (ne_of_gt hsb), mul_one]

theorem continuousAt_thetaHatFn {θi ci cj p : ℝ} (hθi : 0 ≤ θi) (hci : 0 ≤ ci)
    (hthr : ci < θi ^ 4) : ContinuousAt (thetaHatFn ci cj) (rhoSq θi ci, p) := by
  have hb : 0 < betaSq θi ci := betaSq_pos hci hthr
  have hden : betaHat (rhoSq θi ci) ci ≠ 0 := by
    rw [betaHat_rhoSq hθi hci hthr]
    exact ne_of_gt (Real.sqrt_pos.mpr hb)
  have hnum : ContinuousAt (fun sp : ℝ × ℝ => Real.sqrt (sp.2 - cj)) (rhoSq θi ci, p) := by
    have h := ContinuousAt.comp (g := Real.sqrt) (f := fun sp : ℝ × ℝ => sp.2 - cj)
      (x := (rhoSq θi ci, p)) Real.continuous_sqrt.continuousAt
      (continuousAt_snd.sub continuousAt_const)
    simpa [Function.comp_def] using h
  have hde : ContinuousAt (fun sp : ℝ × ℝ => betaHat sp.1 ci) (rhoSq θi ci, p) := by
    have h := ContinuousAt.comp (g := fun s => betaHat s ci)
      (f := fun sp : ℝ × ℝ => sp.1) (x := (rhoSq θi ci, p))
      (continuousAt_betaHat hθi hci hthr) continuousAt_fst
    simpa [Function.comp_def] using h
  exact hnum.div hde hden

/-! ### `‖A a‖²` as a dot product

`Matrix.mulVec` returns a plain `Fin r → ℝ`, whose `‖·‖` is the sup norm of the pi type. The
Euclidean square norm is the dot product. -/

/-- `‖A a‖²`. -/
noncomputable def sqNormMulVec {r D : ℕ} (A : Matrix (Fin r) (Fin D) ℝ)
    (a : EuclideanSpace ℝ (Fin D)) : ℝ :=
  (A *ᵥ WithLp.ofLp a) ⬝ᵥ (A *ᵥ WithLp.ofLp a)

theorem sqNormMulVec_eq_norm_sq {r D : ℕ} (A : Matrix (Fin r) (Fin D) ℝ)
    (a : EuclideanSpace ℝ (Fin D)) :
    sqNormMulVec A a
      = ‖(WithLp.toLp 2 (A *ᵥ WithLp.ofLp a) : EuclideanSpace ℝ (Fin r))‖ ^ 2 := by
  have h := inner_euclidean_eq_dotProduct
    (WithLp.toLp 2 (A *ᵥ WithLp.ofLp a) : EuclideanSpace ℝ (Fin r))
    (WithLp.toLp 2 (A *ᵥ WithLp.ofLp a) : EuclideanSpace ℝ (Fin r))
  rw [real_inner_self_eq_norm_sq] at h
  exact h.symm

/-! ### `topDir`: a measurable unit top right singular vector -/

section TopDir

variable {D : ℕ}

/-- A measurable unit top right singular vector of `X`: the normalized first nonzero column of
`topProj (Xᵀ X)`. It is `delocDir` at `v = 0`, because `perpOf 0` is the identity. -/
noncomputable def topDir {r : ℕ} (X : Matrix (Fin r) (Fin D) ℝ) :
    EuclideanSpace ℝ (Fin D) :=
  delocDir (0 : EuclideanSpace ℝ (Fin D)) X

theorem perpTopVec_zero {r : ℕ} (X : Matrix (Fin r) (Fin D) ℝ) (k : Fin D) :
    perpTopVec (0 : EuclideanSpace ℝ (Fin D)) X k
      = topProj (Xᵀ * X) (isHermitian_transpose_mul_self X)
        (EuclideanSpace.single k (1 : ℝ)) := by
  simp [perpTopVec, perpOf]

theorem measurable_topDir {r : ℕ} :
    Measurable (fun X : Matrix (Fin r) (Fin D) ℝ => topDir X) :=
  measurable_delocDir 0

/-- The top eigenprojector of a nonzero space has a nonzero column. -/
private theorem exists_topProj_single_ne_zero {r : ℕ} (hD : 0 < D)
    (X : Matrix (Fin r) (Fin D) ℝ) :
    ∃ k, topProj (Xᵀ * X) (isHermitian_transpose_mul_self X)
      (EuclideanSpace.single k (1 : ℝ)) ≠ 0 := by
  by_contra hcon
  push Not at hcon
  have hP : ∀ w : EuclideanSpace ℝ (Fin D),
      topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) w
        = (topSpace (Xᵀ * X) (isHermitian_transpose_mul_self X)).starProjection w :=
    fun _ => rfl
  have hzero : ∀ w : EuclideanSpace ℝ (Fin D),
      topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) w = 0 := by
    intro w
    refine PiLp.ext fun k => ?_
    have h1 : (topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) w) k
        = ⟪EuclideanSpace.single k (1 : ℝ),
            topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) w⟫_ℝ := by
      rw [EuclideanSpace.inner_single_left]
      simp
    have h2 : ⟪EuclideanSpace.single k (1 : ℝ),
          topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) w⟫_ℝ
        = ⟪topProj (Xᵀ * X) (isHermitian_transpose_mul_self X)
            (EuclideanSpace.single k (1 : ℝ)), w⟫_ℝ := by
      rw [hP, hP]
      exact (Submodule.inner_starProjection_left_eq_right _ _ _).symm
    rw [h1, h2, hcon k, inner_zero_left]
    simp
  have hmem : topProj (Xᵀ * X) (isHermitian_transpose_mul_self X)
      (vMax (Xᵀ * X) (isHermitian_transpose_mul_self X))
      = vMax (Xᵀ * X) (isHermitian_transpose_mul_self X) :=
    Submodule.starProjection_eq_self_iff.mpr
      (mem_topSpace_vMax hD (Xᵀ * X) (isHermitian_transpose_mul_self X))
  have hn : ‖vMax (Xᵀ * X) (isHermitian_transpose_mul_self X)‖ = 1 :=
    norm_vMax hD (Xᵀ * X) (isHermitian_transpose_mul_self X)
  rw [hzero] at hmem
  rw [← hmem, norm_zero] at hn
  exact zero_ne_one hn

theorem mem_topSpace_topDir {r : ℕ} (X : Matrix (Fin r) (Fin D) ℝ) :
    topDir X ∈ topSpace (Xᵀ * X) (isHermitian_transpose_mul_self X) := by
  rcases delocDir_spec (0 : EuclideanSpace ℝ (Fin D)) X with ⟨h0, -⟩ | ⟨k, -, hdir⟩
  · rw [topDir, h0]
    exact Submodule.zero_mem _
  · rw [topDir, hdir, perpTopVec_zero]
    exact Submodule.smul_mem _ _ (Submodule.starProjection_apply_mem _ _)

theorem norm_topDir {r : ℕ} (hD : 0 < D) (X : Matrix (Fin r) (Fin D) ℝ) : ‖topDir X‖ = 1 := by
  rcases delocDir_spec (0 : EuclideanSpace ℝ (Fin D)) X with ⟨-, hall⟩ | ⟨k, hk, hdir⟩
  · obtain ⟨k, hk⟩ := exists_topProj_single_ne_zero hD X
    exact absurd ((perpTopVec_zero X k).symm.trans (hall k)) hk
  · rw [topDir, hdir, norm_smul, norm_inv, Real.norm_eq_abs, abs_of_nonneg (norm_nonneg _)]
    exact inv_mul_cancel₀ (norm_ne_zero_iff.mpr hk)

theorem dotProduct_topDir_self {r : ℕ} (hD : 0 < D) (X : Matrix (Fin r) (Fin D) ℝ) :
    WithLp.ofLp (topDir X) ⬝ᵥ WithLp.ofLp (topDir X) = 1 :=
  R2.dotProduct_ofLp_self (norm_topDir hD X)

/-- On the simple event, `topDir X` is a unit top eigenvector, so it computes `overlap`. -/
theorem overlap_eq_inner_topDir {r : ℕ} (hD : 0 < D) (X : Matrix (Fin r) (Fin D) ℝ)
    (w : EuclideanSpace ℝ (Fin D))
    (hsimple : TopSimple (Xᵀ * X) (isHermitian_transpose_mul_self X)) :
    overlap X w = ⟪topDir X, w⟫_ℝ ^ 2 :=
  overlap_eq_inner_sq X w hsimple (mem_topSpace_topDir X) (norm_topDir hD X)

end TopDir

/-! ### The law of `Z a` for a unit `a`, and the two second moments -/

section Transport

variable {r D : ℕ}

/-- **The transport step of item P.** For a unit `a`, the column combination `Z ↦ Z a` carries
the canonical Gaussian matrix law to the standard Gaussian vector law of `ℝ^r`. -/
theorem measurePreserving_mulVec_unit (hD : 0 < D) {a : Fin D → ℝ} (ha : a ⬝ᵥ a = 1) :
    MeasurePreserving (fun Z : Matrix (Fin r) (Fin D) ℝ => Z *ᵥ a)
      (gaussianMatrix r D) (R2.piGauss r) := by
  classical
  set j0 : Fin D := ⟨0, hD⟩ with hj0
  set aE : EuclideanSpace ℝ (Fin D) := WithLp.toLp 2 a with haEdef
  have haE : ‖aE‖ = 1 := by
    have h : ⟪aE, aE⟫_ℝ = a ⬝ᵥ a := inner_euclidean_eq_dotProduct aE aE
    rw [real_inner_self_eq_norm_sq, ha] at h
    have h2 : ‖aE‖ = Real.sqrt (‖aE‖ ^ 2) := (Real.sqrt_sq (norm_nonneg _)).symm
    rw [h2, h, Real.sqrt_one]
  have hcard : Module.finrank ℝ (EuclideanSpace ℝ (Fin D)) = Fintype.card (Fin D) := by
    rw [finrank_euclideanSpace_fin, Fintype.card_fin]
  have hortho : Orthonormal ℝ
      (({j0} : Set (Fin D)).domRestrict (fun _ : Fin D => aE)) := by
    refine ⟨fun _ => haE, ?_⟩
    intro x y hxy
    exact absurd (Subtype.ext ((x.2 : (x : Fin D) = j0).trans
      (y.2 : (y : Fin D) = j0).symm)) hxy
  obtain ⟨b, hb⟩ := hortho.exists_orthonormalBasis_extension_of_card_eq hcard
  have hb0 : WithLp.ofLp (b j0) = a := by
    rw [hb j0 rfl]
  set U : Matrix (Fin D) (Fin D) ℝ := Matrix.of fun i i' => WithLp.ofLp (b i) i' with hUdef
  have hUUt : U * Uᵀ = 1 := by
    ext i i'
    have hb' := (orthonormal_iff_ite (𝕜 := ℝ)).1 b.orthonormal i i'
    rw [inner_euclidean_eq_dotProduct] at hb'
    rw [Matrix.one_apply, ← hb']
    simp [hUdef, Matrix.mul_apply, dotProduct]
  have hU : Uᵀ * U = 1 := mul_eq_one_comm.1 hUUt
  have hrow : ∀ l, U j0 l = a l := fun l => congrFun hb0 l
  have hUj0 : U j0 = a := funext hrow
  have h1 := measurePreserving_rowMulVec hU r
  have h2pi : MeasurePreserving
      (fun (Y : Fin r → Fin D → ℝ) (i : Fin r) => Y i j0)
      (Measure.pi fun _ : Fin r => Measure.pi fun _ : Fin D => gaussianReal 0 1)
      (Measure.pi fun _ : Fin r => gaussianReal 0 1) :=
    measurePreserving_pi _ _ fun _ => measurePreserving_eval _ j0
  have h2 : MeasurePreserving
      (fun (Y : Matrix (Fin r) (Fin D) ℝ) (i : Fin r) => Y i j0)
      (gaussianMatrix r D) (R2.piGauss r) := h2pi
  have hcomp := h2.comp h1
  have hfun : (fun Z : Matrix (Fin r) (Fin D) ℝ => Z *ᵥ a)
      = (fun (Y : Matrix (Fin r) (Fin D) ℝ) (i : Fin r) => Y i j0) ∘
        (fun (Y : Matrix (Fin r) (Fin D) ℝ) (i : Fin r) => U *ᵥ Y i) := by
    funext Z i
    change ∑ l, Z i l * a l = ∑ l, U j0 l * Z i l
    exact Finset.sum_congr rfl fun l _ => by rw [hrow l]; ring
  rw [hfun]
  exact hcomp

theorem dotProduct_self_eq_sum_sq (y : Fin r → ℝ) : y ⬝ᵥ y = ∑ k, y k ^ 2 :=
  Finset.sum_congr rfl fun k _ => (sq (y k)).symm

/-- The centered second moment of `∑ y_k²` under the standard Gaussian vector law. -/
theorem integral_sq_sum_sq_sub (r : ℕ) :
    ∫ y, ((∑ k, y k ^ 2) - (r : ℝ)) ^ 2 ∂(R2.piGauss r) = R2.varSq * r := by
  have h := R2.Centered.integral_sq_sum (d := r) R2.centered_sq_sub_one (fun _ => (1 : ℝ))
  have hL : ∀ y : Fin r → ℝ, (∑ k, (1 : ℝ) * ((y k) ^ 2 - 1)) = (∑ k, y k ^ 2) - (r : ℝ) := by
    intro y
    simp [Finset.sum_sub_distrib]
  simp only [hL] at h
  rw [h]
  simp [R2.varSq]

theorem integrable_sq_sum_sq_sub (r : ℕ) :
    Integrable (fun y : Fin r → ℝ => ((∑ k, y k ^ 2) - (r : ℝ)) ^ 2) (R2.piGauss r) := by
  have h := R2.Centered.integrable_sq_sum (d := r) R2.centered_sq_sub_one (fun _ => (1 : ℝ))
  have hL : ∀ y : Fin r → ℝ, (∑ k, (1 : ℝ) * ((y k) ^ 2 - 1)) = (∑ k, y k ^ 2) - (r : ℝ) := by
    intro y
    simp [Finset.sum_sub_distrib]
  simpa only [hL] using h

/-- The pointwise identity behind the two moment computations. -/
private theorem center_eq (D r : ℕ) (y : Fin r → ℝ) :
    ((D : ℝ)⁻¹ * (y ⬝ᵥ y) - (r : ℝ) / D) ^ 2
      = ((D : ℝ)⁻¹) ^ 2 * ((∑ k, y k ^ 2) - (r : ℝ)) ^ 2 := by
  rw [dotProduct_self_eq_sum_sq, div_eq_inv_mul, ← mul_sub, mul_pow]

/-- **Item P, the second moment.** `∫ (d⁻¹ ‖Z a‖² - r/d)² = varSq · r / d²`, uniform over unit
vectors `a`. -/
theorem integral_sq_normSqMulVec (hD : 0 < D) {a : Fin D → ℝ} (ha : a ⬝ᵥ a = 1) :
    ∫ Z, ((D : ℝ)⁻¹ * ((Z *ᵥ a) ⬝ᵥ (Z *ᵥ a)) - (r : ℝ) / D) ^ 2 ∂(gaussianMatrix r D)
      = R2.varSq * r / (D : ℝ) ^ 2 := by
  have hmp := measurePreserving_mulVec_unit (r := r) hD ha
  have hint : Integrable
      (fun y : Fin r → ℝ => ((D : ℝ)⁻¹ * (y ⬝ᵥ y) - (r : ℝ) / D) ^ 2) (R2.piGauss r) := by
    simp only [center_eq]
    exact (integrable_sq_sum_sq_sub r).const_mul _
  have htr : ∫ Z, ((D : ℝ)⁻¹ * ((Z *ᵥ a) ⬝ᵥ (Z *ᵥ a)) - (r : ℝ) / D) ^ 2
        ∂(gaussianMatrix r D)
      = ∫ y, ((D : ℝ)⁻¹ * (y ⬝ᵥ y) - (r : ℝ) / D) ^ 2 ∂(R2.piGauss r) := by
    have h := integral_map (φ := fun Z : Matrix (Fin r) (Fin D) ℝ => Z *ᵥ a)
      (μ := gaussianMatrix r D) hmp.measurable.aemeasurable
      (f := fun y : Fin r → ℝ => ((D : ℝ)⁻¹ * (y ⬝ᵥ y) - (r : ℝ) / D) ^ 2)
      (by rw [hmp.map_eq]; exact hint.aestronglyMeasurable)
    rw [hmp.map_eq] at h
    exact h.symm
  rw [htr]
  simp only [center_eq]
  rw [integral_const_mul, integral_sq_sum_sq_sub, inv_pow]
  ring

theorem integrable_sq_normSqMulVec (hD : 0 < D) {a : Fin D → ℝ} (ha : a ⬝ᵥ a = 1) :
    Integrable
      (fun Z : Matrix (Fin r) (Fin D) ℝ =>
        ((D : ℝ)⁻¹ * ((Z *ᵥ a) ⬝ᵥ (Z *ᵥ a)) - (r : ℝ) / D) ^ 2) (gaussianMatrix r D) := by
  have hmp := measurePreserving_mulVec_unit (r := r) hD ha
  have hint : Integrable
      (fun y : Fin r → ℝ => ((D : ℝ)⁻¹ * (y ⬝ᵥ y) - (r : ℝ) / D) ^ 2) (R2.piGauss r) := by
    simp only [center_eq]
    exact (integrable_sq_sum_sq_sub r).const_mul _
  exact memLp_one_iff_integrable.mp
    ((memLp_one_iff_integrable.mpr hint).comp_measurePreserving hmp)

/-- **The cross term, second moment.** `∫ (x ⬝ᵥ Z a)² = 1` for unit `x` and unit `a`. -/
theorem integral_sq_dotProduct_mulVec (hD : 0 < D) {x : Fin r → ℝ} (hx : x ⬝ᵥ x = 1)
    {a : Fin D → ℝ} (ha : a ⬝ᵥ a = 1) :
    ∫ Z, (x ⬝ᵥ (Z *ᵥ a)) ^ 2 ∂(gaussianMatrix r D) = 1 := by
  have hmp := measurePreserving_mulVec_unit (r := r) hD ha
  have hsum : ∀ y : Fin r → ℝ, x ⬝ᵥ y = ∑ k, x k * y k := fun _ => rfl
  have hint : Integrable (fun y : Fin r → ℝ => (x ⬝ᵥ y) ^ 2) (R2.piGauss r) := by
    simp only [hsum]
    exact R2.Centered.integrable_sq_sum R2.centered_id x
  have htr : ∫ Z, (x ⬝ᵥ (Z *ᵥ a)) ^ 2 ∂(gaussianMatrix r D)
      = ∫ y, (x ⬝ᵥ y) ^ 2 ∂(R2.piGauss r) := by
    have h := integral_map (φ := fun Z : Matrix (Fin r) (Fin D) ℝ => Z *ᵥ a)
      (μ := gaussianMatrix r D) hmp.measurable.aemeasurable
      (f := fun y : Fin r → ℝ => (x ⬝ᵥ y) ^ 2)
      (by rw [hmp.map_eq]; exact hint.aestronglyMeasurable)
    rw [hmp.map_eq] at h
    exact h.symm
  rw [htr]
  simp only [hsum]
  rw [R2.Centered.integral_sq_sum R2.centered_id x, R2.integral_sq_gauss]
  have hx2 : ∑ k, x k ^ 2 = x ⬝ᵥ x := (dotProduct_self_eq_sum_sq x).symm
  rw [hx2, hx, one_mul]

theorem integrable_sq_dotProduct_mulVec (hD : 0 < D) {x : Fin r → ℝ}
    {a : Fin D → ℝ} (ha : a ⬝ᵥ a = 1) :
    Integrable (fun Z : Matrix (Fin r) (Fin D) ℝ => (x ⬝ᵥ (Z *ᵥ a)) ^ 2)
      (gaussianMatrix r D) := by
  have hmp := measurePreserving_mulVec_unit (r := r) hD ha
  have hsum : ∀ y : Fin r → ℝ, x ⬝ᵥ y = ∑ k, x k * y k := fun _ => rfl
  have hint : Integrable (fun y : Fin r → ℝ => (x ⬝ᵥ y) ^ 2) (R2.piGauss r) := by
    simp only [hsum]
    exact R2.Centered.integrable_sq_sum R2.centered_id x
  exact memLp_one_iff_integrable.mp
    ((memLp_one_iff_integrable.mpr hint).comp_measurePreserving hmp)

/-! ### Measurability of the two kernels -/

theorem measurable_mulVec_pair :
    Measurable fun q : (Fin D → ℝ) × Matrix (Fin r) (Fin D) ℝ => q.2 *ᵥ q.1 := by
  refine measurable_pi_lambda _ fun k => ?_
  simp only [Matrix.mulVec, dotProduct]
  refine Finset.measurable_sum _ fun l _ => ?_
  exact ((measurable_pi_apply l).comp
    ((measurable_pi_apply k).comp measurable_snd)).mul
      ((measurable_pi_apply l).comp measurable_fst)

theorem measurable_dotProduct_pair :
    Measurable fun q : (Fin D → ℝ) × Matrix (Fin r) (Fin D) ℝ =>
      (q.2 *ᵥ q.1) ⬝ᵥ (q.2 *ᵥ q.1) := by
  simp only [dotProduct]
  refine Finset.measurable_sum _ fun k _ => ?_
  exact ((measurable_pi_apply k).comp measurable_mulVec_pair).mul
    ((measurable_pi_apply k).comp measurable_mulVec_pair)

theorem measurable_dotProduct_left_pair (x : Fin r → ℝ) :
    Measurable fun q : (Fin D → ℝ) × Matrix (Fin r) (Fin D) ℝ => x ⬝ᵥ (q.2 *ᵥ q.1) := by
  simp only [dotProduct]
  refine Finset.measurable_sum _ fun k _ => ?_
  exact measurable_const.mul ((measurable_pi_apply k).comp measurable_mulVec_pair)

end Transport

end ThetaEst

/-! ### Item P for one table -/

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- `‖E a‖² = d⁻¹ ‖Z a‖²`. -/
theorem sqNormMulVec_E (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N)
    (a : EuclideanSpace ℝ (Fin (d N))) :
    ThetaEst.sqNormMulVec (m.E N ω) a
      = ((d N : ℝ))⁻¹ *
        ((m.Z N ω *ᵥ WithLp.ofLp a) ⬝ᵥ (m.Z N ω *ᵥ WithLp.ofLp a)) := by
  change (((Real.sqrt (d N))⁻¹ • m.Z N ω) *ᵥ WithLp.ofLp a) ⬝ᵥ
      (((Real.sqrt (d N))⁻¹ • m.Z N ω) *ᵥ WithLp.ofLp a) = _
  rw [Matrix.smul_mulVec, smul_dotProduct, dotProduct_smul, smul_eq_mul, smul_eq_mul,
    ← mul_assoc, R2.sqrt_inv_mul_sqrt_inv]

/-- The paper's expansion of `‖X a‖²` (`main_paper.tex:2470`). -/
theorem sqNormMulVec_X (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N)
    (a : EuclideanSpace ℝ (Fin (d N))) :
    ThetaEst.sqNormMulVec (m.X N ω) a
      = m.θ ^ 2 * ⟪m.v N, a⟫_ℝ ^ 2
        + 2 * (m.θ * ⟪m.v N, a⟫_ℝ) *
          (WithLp.ofLp (m.u N) ⬝ᵥ (m.E N ω *ᵥ WithLp.ofLp a))
        + ThetaEst.sqNormMulVec (m.E N ω) a := by
  have hinner : ⟪m.v N, a⟫_ℝ = WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp a :=
    inner_euclidean_eq_dotProduct (m.v N) a
  have hx : WithLp.ofLp (m.u N) ⬝ᵥ WithLp.ofLp (m.u N) = 1 := m.dotProduct_u_self N
  have hvm : Matrix.vecMulVec (WithLp.ofLp (m.u N)) (WithLp.ofLp (m.v N)) *ᵥ WithLp.ofLp a
      = (WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp a) • WithLp.ofLp (m.u N) := by
    funext k
    change ∑ l, Matrix.vecMulVec (WithLp.ofLp (m.u N)) (WithLp.ofLp (m.v N)) k l
        * WithLp.ofLp a l
      = (WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp a) * WithLp.ofLp (m.u N) k
    rw [dotProduct, Finset.sum_mul]
    refine Finset.sum_congr rfl fun l _ => ?_
    rw [Matrix.vecMulVec_apply]
    ring
  have hX : m.X N ω *ᵥ WithLp.ofLp a
      = (m.θ * (WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp a)) • WithLp.ofLp (m.u N)
        + m.E N ω *ᵥ WithLp.ofLp a := by
    change (m.θ • Matrix.vecMulVec (WithLp.ofLp (m.u N)) (WithLp.ofLp (m.v N)) + m.E N ω)
        *ᵥ WithLp.ofLp a = _
    rw [Matrix.add_mulVec, Matrix.smul_mulVec, hvm, smul_smul]
  have hexp : ∀ (s : ℝ) (x y : Fin (n N) → ℝ), x ⬝ᵥ x = 1 →
      (s • x + y) ⬝ᵥ (s • x + y) = s ^ 2 + 2 * s * (x ⬝ᵥ y) + y ⬝ᵥ y := by
    intro s x y hxx
    simp only [add_dotProduct, dotProduct_add, smul_dotProduct, dotProduct_smul, smul_eq_mul,
      dotProduct_comm y x, hxx]
    ring
  change (m.X N ω *ᵥ WithLp.ofLp a) ⬝ᵥ (m.X N ω *ᵥ WithLp.ofLp a) = _
  rw [hX, hexp _ _ _ hx, hinner]
  change _ = m.θ ^ 2 * (WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp a) ^ 2
      + 2 * (m.θ * (WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp a)) *
        (WithLp.ofLp (m.u N) ⬝ᵥ (m.E N ω *ᵥ WithLp.ofLp a))
      + (m.E N ω *ᵥ WithLp.ofLp a) ⬝ᵥ (m.E N ω *ᵥ WithLp.ofLp a)
  ring

/-- **Item P**, `lem:noise_projection_concentration` for Gaussian noise and a deterministic
direction: `‖E_N a_N‖² → c` in probability. -/
theorem noise_projection_tendsto (m : SpikedModel μ n d) {c : ℝ}
    (hG : m.GaussianNoise) (hreg : m.Regime c)
    (a : (N : ℕ) → EuclideanSpace ℝ (Fin (d N))) (ha : ∀ N, ‖a N‖ = 1) :
    TendstoInProb μ (fun N ω => ThetaEst.sqNormMulVec (m.E N ω) (a N)) c := by
  have hd : Tendsto d atTop atTop := hreg.2.1
  have hbnd : ∀ ε : ℝ, 0 < ε → ∀ N,
      μ N {ω | ε ≤ |ThetaEst.sqNormMulVec (m.E N ω) (a N) - (n N : ℝ) / d N|}
        ≤ ENNReal.ofReal (R2.varSq * (n N) / (d N : ℝ) ^ 2 / ε ^ 2) := by
    intro ε hε N
    have hD : 0 < d N := m.hd N
    have haN : WithLp.ofLp (a N) ⬝ᵥ WithLp.ofLp (a N) = 1 := R2.dotProduct_ofLp_self (ha N)
    set G : Matrix (Fin (n N)) (Fin (d N)) ℝ → ℝ := fun Z =>
      |((d N : ℝ))⁻¹ * ((Z *ᵥ WithLp.ofLp (a N)) ⬝ᵥ (Z *ᵥ WithLp.ofLp (a N)))
        - (n N : ℝ) / d N| with hGdef
    have hmeas : Measurable G := by
      rw [hGdef]
      exact (((ThetaEst.measurable_dotProduct_pair.comp
        (measurable_const.prodMk measurable_id)).const_mul _).sub measurable_const).abs
    have hset : {ω | ε ≤ |ThetaEst.sqNormMulVec (m.E N ω) (a N) - (n N : ℝ) / d N|}
        = {ω | ε ≤ G (m.Z N ω)} := by
      ext ω
      rw [Set.mem_ofPred_eq, Set.mem_ofPred_eq, hGdef]
      simp only
      rw [m.sqNormMulVec_E N ω (a N)]
    have key : μ N {ω | ε ≤ G (m.Z N ω)}
        = gaussianMatrix (n N) (d N) {Z | ε ≤ G Z} :=
      (hG N).measure_eq (measurableSet_le measurable_const hmeas)
    rw [hset, key]
    refine R2.meas_ge_le_of_integral_sq _ hmeas ?_ hε ?_
    · rw [hGdef]
      simpa [sq_abs] using ThetaEst.integrable_sq_normSqMulVec (r := n N) hD haN
    · rw [hGdef]
      simpa [sq_abs] using
        le_of_eq (ThetaEst.integral_sq_normSqMulVec (r := n N) hD haN)
  have hK : Tendsto (fun N => R2.varSq * (n N : ℝ) / (d N : ℝ) ^ 2) atTop (𝓝 0) := by
    have hdR : Tendsto (fun N => ((d N : ℝ))) atTop atTop :=
      tendsto_natCast_atTop_atTop.comp hd
    have hinv : Tendsto (fun N => ((d N : ℝ))⁻¹) atTop (𝓝 0) := hdR.inv_tendsto_atTop
    have hprod := (hreg.2.2.const_mul R2.varSq).mul hinv
    have hz : R2.varSq * c * 0 = 0 := by ring
    rw [hz] at hprod
    exact hprod.congr fun N => by ring
  have h0 := ThetaEst.tendstoInProb_of_meas_le hK hbnd
  have hfin := h0.add_tendsto hreg.2.2
  simpa using hfin

end SpikedModel

/-! ### The estimator and `thm:theta_est` -/

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

/-- The two cross-table concentration limits that `thm:theta_est` consumes: the noise of
table `j` against the top direction of table `i`. Both hold for independent Gaussian tables
(`thetaEstLaw_of_gaussian`); they are the only Gaussian input of `thm_theta_est` (F31).
* `noiseProj`: item P with a random direction, `‖E_j v̂_i‖² → c_j`.
* `crossTerm`: `u_jᵀ E_j v̂_i → 0`. -/
structure ThetaEstLaw (m : MultiTableModel μ M n d) (i j : Fin M) (cj : ℝ) : Prop where
  noiseProj : TendstoInProb μ (fun N ω => ThetaEst.sqNormMulVec ((m.tbl j).E N ω)
    (ThetaEst.topDir ((m.tbl i).X N ω))) cj
  crossTerm : TendstoInProb μ (fun N ω => WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ
    ((m.tbl j).E N ω *ᵥ WithLp.ofLp (ThetaEst.topDir ((m.tbl i).X N ω)))) 0

/-- **The Fubini step.** A second moment bound in the noise of table `j`, uniform over unit
directions, and a direction that is a measurable function of the noise of table `i`, give a
Chebyshev bound for the random direction. Independence across tables enters through `hG`; the
route is that of `measure_deloc_le`, with the two indices swapped. -/
theorem measure_randomDir_le (m : MultiTableModel μ M n d) (hG : m.JointGaussianNoise)
    {i j : Fin M} (hij : i ≠ j) (N : ℕ)
    (A : Matrix (Fin (n i N)) (Fin (d N)) ℝ → (Fin (d N) → ℝ)) (hAmeas : Measurable A)
    (hAunit : ∀ Z, A Z ⬝ᵥ A Z = 1)
    (F : (Fin (d N) → ℝ) → Matrix (Fin (n j N)) (Fin (d N)) ℝ → ℝ)
    (hFmeas : Measurable fun q : (Fin (d N) → ℝ) × Matrix (Fin (n j N)) (Fin (d N)) ℝ =>
      F q.1 q.2)
    (hFint : ∀ a : Fin (d N) → ℝ, a ⬝ᵥ a = 1 →
      Integrable (fun Z => F a Z ^ 2) (gaussianMatrix (n j N) (d N)))
    {C ε : ℝ} (hε : 0 < ε)
    (hC : ∀ a : Fin (d N) → ℝ, a ⬝ᵥ a = 1 →
      ∫ Z, F a Z ^ 2 ∂(gaussianMatrix (n j N) (d N)) ≤ C) :
    μ N {ω | ε ≤ |F (A ((m.tbl i).Z N ω)) ((m.tbl j).Z N ω)|}
      ≤ ENNReal.ofReal (C / ε ^ 2) := by
  classical
  have hprob : ∀ k : Fin M, IsProbabilityMeasure (gaussianMatrix (n k N) (d N)) := fun k =>
    isProbabilityMeasure_gaussianMatrix _ _
  have hinst := hprob
  have hSmeas : MeasurableSet {q : Matrix (Fin (n j N)) (Fin (d N)) ℝ ×
      Matrix (Fin (n i N)) (Fin (d N)) ℝ | ε ≤ |F (A q.2) q.1|} := by
    have h1 : Measurable fun q : Matrix (Fin (n j N)) (Fin (d N)) ℝ ×
        Matrix (Fin (n i N)) (Fin (d N)) ℝ => F (A q.2) q.1 :=
      hFmeas.comp ((hAmeas.comp measurable_snd).prodMk measurable_fst)
    exact measurableSet_le measurable_const h1.abs
  have hpairmeas : MeasurableSet {Zs : (k : Fin M) → Matrix (Fin (n k N)) (Fin (d N)) ℝ |
      (Zs j, Zs i) ∈ {q : Matrix (Fin (n j N)) (Fin (d N)) ℝ ×
        Matrix (Fin (n i N)) (Fin (d N)) ℝ | ε ≤ |F (A q.2) q.1|}} :=
    ((measurable_pi_apply j).prodMk (measurable_pi_apply i)) hSmeas
  have key : μ N {ω | ε ≤ |F (A ((m.tbl i).Z N ω)) ((m.tbl j).Z N ω)|}
      = (Measure.pi fun k : Fin M => gaussianMatrix (n k N) (d N))
        {Zs | (Zs j, Zs i) ∈ {q : Matrix (Fin (n j N)) (Fin (d N)) ℝ ×
          Matrix (Fin (n i N)) (Fin (d N)) ℝ | ε ≤ |F (A q.2) q.1|}} :=
    (hG N).measure_eq hpairmeas
  rw [key]
  refine pi_pair_le (fun k : Fin M => gaussianMatrix (n k N) (d N)) hprob hij.symm hSmeas _ ?_
  intro zi
  have hmeas' : Measurable fun zj : Matrix (Fin (n j N)) (Fin (d N)) ℝ => |F (A zi) zj| :=
    (hFmeas.comp (measurable_const.prodMk measurable_id)).abs
  have hint' : Integrable (fun zj => |F (A zi) zj| ^ 2) (gaussianMatrix (n j N) (d N)) := by
    simpa [sq_abs] using hFint (A zi) (hAunit zi)
  have hC' : ∫ zj, |F (A zi) zj| ^ 2 ∂(gaussianMatrix (n j N) (d N)) ≤ C := by
    simpa [sq_abs] using hC (A zi) (hAunit zi)
  exact R2.meas_ge_le_of_integral_sq _ hmeas' hint' hε hC'

/-- Item P with the direction read off table `i`: `‖E_j v̂_i‖² → c_j`. -/
theorem noise_projection_topDir_tendsto (m : MultiTableModel μ M n d)
    (hG : m.JointGaussianNoise) {i j : Fin M} (hij : i ≠ j) {cj : ℝ}
    (hregj : (m.tbl j).Regime cj) :
    TendstoInProb μ (fun N ω => ThetaEst.sqNormMulVec ((m.tbl j).E N ω)
      (ThetaEst.topDir ((m.tbl i).X N ω))) cj := by
  have hd : Tendsto d atTop atTop := hregj.2.1
  have hbnd : ∀ ε : ℝ, 0 < ε → ∀ N,
      μ N {ω | ε ≤ |ThetaEst.sqNormMulVec ((m.tbl j).E N ω)
          (ThetaEst.topDir ((m.tbl i).X N ω)) - (n j N : ℝ) / d N|}
        ≤ ENNReal.ofReal (R2.varSq * (n j N) / (d N : ℝ) ^ 2 / ε ^ 2) := by
    intro ε hε N
    have hD : 0 < d N := (m.tbl j).hd N
    set A : Matrix (Fin (n i N)) (Fin (d N)) ℝ → (Fin (d N) → ℝ) :=
      fun Z => WithLp.ofLp (ThetaEst.topDir (m.dataOf i N Z)) with hAdef
    set F : (Fin (d N) → ℝ) → Matrix (Fin (n j N)) (Fin (d N)) ℝ → ℝ :=
      fun a Z => ((d N : ℝ))⁻¹ * ((Z *ᵥ a) ⬝ᵥ (Z *ᵥ a)) - (n j N : ℝ) / d N with hFdef
    have hAmeas : Measurable A := by
      rw [hAdef]
      exact (WithLp.measurable_ofLp 2 (Fin (d N) → ℝ)).comp
        (ThetaEst.measurable_topDir.comp (m.measurable_dataOf i N))
    have hAunit : ∀ Z, A Z ⬝ᵥ A Z = 1 := by
      intro Z
      rw [hAdef]
      exact ThetaEst.dotProduct_topDir_self hD _
    have hFmeas : Measurable fun q : (Fin (d N) → ℝ) ×
        Matrix (Fin (n j N)) (Fin (d N)) ℝ => F q.1 q.2 := by
      rw [hFdef]
      exact (ThetaEst.measurable_dotProduct_pair.const_mul _).sub measurable_const
    have hFint : ∀ a : Fin (d N) → ℝ, a ⬝ᵥ a = 1 →
        Integrable (fun Z => F a Z ^ 2) (gaussianMatrix (n j N) (d N)) := by
      intro a hau
      rw [hFdef]
      exact ThetaEst.integrable_sq_normSqMulVec (r := n j N) hD hau
    have hset : {ω | ε ≤ |ThetaEst.sqNormMulVec ((m.tbl j).E N ω)
          (ThetaEst.topDir ((m.tbl i).X N ω)) - (n j N : ℝ) / d N|}
        = {ω | ε ≤ |F (A ((m.tbl i).Z N ω)) ((m.tbl j).Z N ω)|} := by
      ext ω
      rw [Set.mem_ofPred_eq, Set.mem_ofPred_eq, hFdef, hAdef]
      simp only
      rw [(m.tbl j).sqNormMulVec_E N ω (ThetaEst.topDir ((m.tbl i).X N ω)),
        m.dataOf_eq i N ω]
    rw [hset]
    refine m.measure_randomDir_le hG hij N A hAmeas hAunit F hFmeas hFint hε ?_
    intro a hau
    rw [hFdef]
    exact le_of_eq (ThetaEst.integral_sq_normSqMulVec (r := n j N) hD hau)
  have hK : Tendsto (fun N => R2.varSq * (n j N : ℝ) / (d N : ℝ) ^ 2) atTop (𝓝 0) := by
    have hdR : Tendsto (fun N => ((d N : ℝ))) atTop atTop :=
      tendsto_natCast_atTop_atTop.comp hd
    have hinv : Tendsto (fun N => ((d N : ℝ))⁻¹) atTop (𝓝 0) := hdR.inv_tendsto_atTop
    have hprod := (hregj.2.2.const_mul R2.varSq).mul hinv
    have hz : R2.varSq * cj * 0 = 0 := by ring
    rw [hz] at hprod
    exact hprod.congr fun N => by ring
  have h0 := ThetaEst.tendstoInProb_of_meas_le hK hbnd
  have hfin := h0.add_tendsto hregj.2.2
  simpa using hfin

/-- The cross term `u_jᵀ E_j v̂_i → 0`. -/
theorem cross_term_tendsto (m : MultiTableModel μ M n d) (hG : m.JointGaussianNoise)
    {i j : Fin M} (hij : i ≠ j) (hd : Tendsto d atTop atTop) :
    TendstoInProb μ (fun N ω => WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ
      ((m.tbl j).E N ω *ᵥ WithLp.ofLp (ThetaEst.topDir ((m.tbl i).X N ω)))) 0 := by
  have hbnd : ∀ ε : ℝ, 0 < ε → ∀ N,
      μ N {ω | ε ≤ |WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ
          ((m.tbl j).E N ω *ᵥ WithLp.ofLp (ThetaEst.topDir ((m.tbl i).X N ω)))
          - (fun _ : ℕ => (0 : ℝ)) N|}
        ≤ ENNReal.ofReal (((d N : ℝ))⁻¹ / ε ^ 2) := by
    intro ε hε N
    have hD : 0 < d N := (m.tbl j).hd N
    have hxu : WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ WithLp.ofLp ((m.tbl j).u N) = 1 :=
      (m.tbl j).dotProduct_u_self N
    have hs2 : ((Real.sqrt (d N))⁻¹) ^ 2 = ((d N : ℝ))⁻¹ := by
      rw [sq]
      exact R2.sqrt_inv_mul_sqrt_inv (d N)
    set A : Matrix (Fin (n i N)) (Fin (d N)) ℝ → (Fin (d N) → ℝ) :=
      fun Z => WithLp.ofLp (ThetaEst.topDir (m.dataOf i N Z)) with hAdef
    set F : (Fin (d N) → ℝ) → Matrix (Fin (n j N)) (Fin (d N)) ℝ → ℝ :=
      fun a Z => (Real.sqrt (d N))⁻¹ * (WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ (Z *ᵥ a)) with hFdef
    have hsq : ∀ (a : Fin (d N) → ℝ) (Z : Matrix (Fin (n j N)) (Fin (d N)) ℝ),
        F a Z ^ 2
          = ((d N : ℝ))⁻¹ * (WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ (Z *ᵥ a)) ^ 2 := by
      intro a Z
      rw [hFdef]
      simp only
      rw [mul_pow, hs2]
    have hAmeas : Measurable A := by
      rw [hAdef]
      exact (WithLp.measurable_ofLp 2 (Fin (d N) → ℝ)).comp
        (ThetaEst.measurable_topDir.comp (m.measurable_dataOf i N))
    have hAunit : ∀ Z, A Z ⬝ᵥ A Z = 1 := by
      intro Z
      rw [hAdef]
      exact ThetaEst.dotProduct_topDir_self hD _
    have hFmeas : Measurable fun q : (Fin (d N) → ℝ) ×
        Matrix (Fin (n j N)) (Fin (d N)) ℝ => F q.1 q.2 := by
      rw [hFdef]
      exact (ThetaEst.measurable_dotProduct_left_pair
        (WithLp.ofLp ((m.tbl j).u N))).const_mul _
    have hFint : ∀ a : Fin (d N) → ℝ, a ⬝ᵥ a = 1 →
        Integrable (fun Z => F a Z ^ 2) (gaussianMatrix (n j N) (d N)) := by
      intro a hau
      simp only [hsq]
      exact (ThetaEst.integrable_sq_dotProduct_mulVec (r := n j N) hD
        (x := WithLp.ofLp ((m.tbl j).u N)) hau).const_mul _
    have hset : {ω | ε ≤ |WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ
          ((m.tbl j).E N ω *ᵥ WithLp.ofLp (ThetaEst.topDir ((m.tbl i).X N ω)))
          - (fun _ : ℕ => (0 : ℝ)) N|}
        = {ω | ε ≤ |F (A ((m.tbl i).Z N ω)) ((m.tbl j).Z N ω)|} := by
      ext ω
      rw [Set.mem_ofPred_eq, Set.mem_ofPred_eq, sub_zero, hFdef, hAdef]
      simp only
      have hE : (m.tbl j).E N ω = (Real.sqrt (d N))⁻¹ • (m.tbl j).Z N ω := rfl
      rw [hE, Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul, m.dataOf_eq i N ω]
    rw [hset]
    refine m.measure_randomDir_le hG hij N A hAmeas hAunit F hFmeas hFint hε ?_
    intro a hau
    simp only [hsq]
    rw [integral_const_mul,
      ThetaEst.integral_sq_dotProduct_mulVec (r := n j N) hD hxu hau, mul_one]
  have hK : Tendsto (fun N => ((d N : ℝ))⁻¹) atTop (𝓝 0) := by
    have hdR : Tendsto (fun N => ((d N : ℝ))) atTop atTop :=
      tendsto_natCast_atTop_atTop.comp hd
    exact hdR.inv_tendsto_atTop
  have h0 := ThetaEst.tendstoInProb_of_meas_le (ctr := fun _ => (0 : ℝ)) hK hbnd
  simpa using h0

/-- The Gaussian discharge of `ThetaEstLaw`: independent Gaussian tables `i ≠ j` in the
regime of table `j` give both limits. -/
theorem thetaEstLaw_of_gaussian (m : MultiTableModel μ M n d) {i j : Fin M} (hij : i ≠ j)
    (hG : m.JointGaussianNoise) {cj : ℝ} (hregj : (m.tbl j).Regime cj) :
    m.ThetaEstLaw i j cj :=
  ⟨m.noise_projection_topDir_tendsto hG hij hregj, m.cross_term_tendsto hG hij hregj.2.1⟩

/-- `eq:theta_estimation`: `θ̂_j` from the top singular value of table `i` and the projection of
table `j` on `v̂_i`. Observable: it reads only `X_i`, `X_j`, `c_i` and `c_j`. -/
noncomputable def thetaHat (m : MultiTableModel μ M n d) (i j : Fin M) (ci cj : ℝ) (N : ℕ)
    (ω : Ω N) : ℝ :=
  ThetaEst.thetaHatFn ci cj (gramLamMax ((m.tbl i).X N ω),
    ThetaEst.sqNormMulVec ((m.tbl j).X N ω) (ThetaEst.topDir ((m.tbl i).X N ω)))

/-- The limit of the projection: `‖X_j v̂_i‖² → θ_j² β_i² + c_j` (`main_paper.tex:2470`). -/
theorem sqNormMulVec_X_tendsto (m : MultiTableModel μ M n d) {i j : Fin M} {ci cj : ℝ}
    (law : (m.tbl i).SingleTableLaw ci) (est : m.ThetaEstLaw i j cj) :
    TendstoInProb μ (fun N ω => ThetaEst.sqNormMulVec ((m.tbl j).X N ω)
        (ThetaEst.topDir ((m.tbl i).X N ω)))
      ((m.tbl j).θ ^ 2 * betaSq (m.tbl i).θ ci + cj) := by
  have hov : TendstoInProb μ
      (fun N ω => ⟪(m.tbl i).v N, ThetaEst.topDir ((m.tbl i).X N ω)⟫_ℝ ^ 2)
      (betaSq (m.tbl i).θ ci) := by
    refine law.align.congr fun N => ?_
    filter_upwards [law.topSimple N] with ω hω
    rw [ThetaEst.overlap_eq_inner_topDir ((m.tbl i).hd N) _ _ hω, real_inner_comm]
  have hA : TendstoInProb μ
      (fun N ω => (m.tbl j).θ ^ 2 *
        ⟪(m.tbl i).v N, ThetaEst.topDir ((m.tbl i).X N ω)⟫_ℝ ^ 2)
      ((m.tbl j).θ ^ 2 * betaSq (m.tbl i).θ ci) :=
    hov.const_mul _
  have hT := est.crossTerm
  have hTabs : TendstoInProb μ (fun N ω => |WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ
      ((m.tbl j).E N ω *ᵥ WithLp.ofLp (ThetaEst.topDir ((m.tbl i).X N ω)))|) 0 := by
    have h := hT.comp_continuous (φ := fun t : ℝ => |t|) continuous_abs.continuousAt
    simpa using h
  have hcross : TendstoInProb μ
      (fun N ω => 2 * ((m.tbl j).θ *
          ⟪(m.tbl i).v N, ThetaEst.topDir ((m.tbl i).X N ω)⟫_ℝ) *
        (WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ
          ((m.tbl j).E N ω *ᵥ WithLp.ofLp (ThetaEst.topDir ((m.tbl i).X N ω))))) 0 := by
    refine TendstoInProb.of_le (g := fun N ω => 2 * (m.tbl j).θ *
      |WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ
        ((m.tbl j).E N ω *ᵥ WithLp.ofLp (ThetaEst.topDir ((m.tbl i).X N ω)))|)
      (fun N => ?_) ?_
    · filter_upwards with ω
      have hb : |⟪(m.tbl i).v N, ThetaEst.topDir ((m.tbl i).X N ω)⟫_ℝ| ≤ 1 := by
        have h := abs_real_inner_le_norm ((m.tbl i).v N)
          (ThetaEst.topDir ((m.tbl i).X N ω))
        rwa [(m.tbl i).hv N, ThetaEst.norm_topDir ((m.tbl i).hd N), mul_one] at h
      have hθ : 0 ≤ (m.tbl j).θ := (m.tbl j).hθ
      have h2 : (0 : ℝ) ≤ |WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ
          ((m.tbl j).E N ω *ᵥ WithLp.ofLp (ThetaEst.topDir ((m.tbl i).X N ω)))| :=
        abs_nonneg _
      rw [sub_zero, abs_mul, abs_mul, abs_two, abs_mul, abs_of_nonneg hθ]
      have h3 : (m.tbl j).θ * |⟪(m.tbl i).v N,
          ThetaEst.topDir ((m.tbl i).X N ω)⟫_ℝ| ≤ (m.tbl j).θ := by
        nlinarith [abs_nonneg (⟪(m.tbl i).v N, ThetaEst.topDir ((m.tbl i).X N ω)⟫_ℝ)]
      nlinarith
    · have h := hTabs.const_mul (2 * (m.tbl j).θ)
      simpa using h
  have hR := est.noiseProj
  have hsum := (hA.add hcross).add hR
  have hval : (m.tbl j).θ ^ 2 * betaSq (m.tbl i).θ ci + 0 + cj
      = (m.tbl j).θ ^ 2 * betaSq (m.tbl i).θ ci + cj := by ring
  rw [hval] at hsum
  refine hsum.congr fun N => ?_
  filter_upwards with ω
  have hvv : (m.tbl j).v N = (m.tbl i).v N := m.hv j i N
  rw [(m.tbl j).sqNormMulVec_X N ω (ThetaEst.topDir ((m.tbl i).X N ω)), hvv]

/-- **`thm:theta_est`** (`main_paper.tex:958`). The estimator of `eq:theta_estimation` is
consistent for `θ_j`, whatever the size of `θ_j`, as soon as table `i` is above its own
detection threshold. -/
theorem thm_theta_est (m : MultiTableModel μ M n d) {i j : Fin M} {ci cj : ℝ}
    (hci : 0 ≤ ci) (hthr : ci < (m.tbl i).θ ^ 4) (law : (m.tbl i).SingleTableLaw ci)
    (est : m.ThetaEstLaw i j cj) :
    TendstoInProb μ (fun N ω => m.thetaHat i j ci cj N ω) (m.tbl j).θ := by
  have hp := m.sqNormMulVec_X_tendsto law est
  have hcont := TendstoInProb.comp_continuous₂ (φ := ThetaEst.thetaHatFn ci cj)
    (ThetaEst.continuousAt_thetaHatFn (θi := (m.tbl i).θ) (cj := cj)
      (p := (m.tbl j).θ ^ 2 * betaSq (m.tbl i).θ ci + cj) (m.tbl i).hθ hci hthr)
    law.lamMax hp
  rwa [ThetaEst.thetaHatFn_limit (m.tbl i).hθ (m.tbl j).hθ hci hthr] at hcont

/-- The Gaussian corollary: `SingleTableLaw` of table `i` is discharged by
`singleTableLaw_of_gaussian` and `ThetaEstLaw` by `thetaEstLaw_of_gaussian`. -/
theorem thm_theta_est_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) {i j : Fin M} (hij : i ≠ j) {ci cj : ℝ} (hci : 0 < ci)
    (hthr : ci < (m.tbl i).θ ^ 4) (hG : m.JointGaussianNoise)
    (hregi : (m.tbl i).Regime ci) (hregj : (m.tbl j).Regime cj) :
    TendstoInProb μ (fun N ω => m.thetaHat i j ci cj N ω) (m.tbl j).θ :=
  m.thm_theta_est hci.le hthr
    (SpikedModel.singleTableLaw_of_gaussian hci (m.tbl i) hregi (m.gaussianNoise_of_joint hG i))
    (m.thetaEstLaw_of_gaussian hij hG hregj)

/-! ### The general-noise discharge of `ThetaEstLaw` (Stage 0 of `notes/NONGAUSSIAN_SCOPE.md`)

The same two limits for a fixed i.i.d. noise law `ν` with mean `0`, variance `1` and a
finite fourth moment (`NoiseLaw ν`, `Prob/NoiseLaw.lean`); `assum:general_noise` of the paper
(`main_paper.tex:244`). The Fubini step is stated for a family of table laws so that both the
Gaussian and the general case instantiate it; the moment inputs are those of
`Prob/NoiseMoments.lean`, with the bound `(ν₄ + 2) n_j / d²` in place of the Gaussian identity
`varSq · n_j / d²`. Nothing above this line changes; `thetaEstLaw_of_gaussian` keeps its own
proof, and the `example` below records that the general theorem gives it back. -/

/-- **The Fubini step, general form.** `measure_randomDir_le` for an arbitrary family of
table laws `L`, with the joint law `Measure.pi L` at the fixed `N`. -/
theorem measure_randomDir_le_general (m : MultiTableModel μ M n d) {i j : Fin M} (hij : i ≠ j)
    (N : ℕ) (L : ∀ k : Fin M, Measure (Matrix (Fin (n k N)) (Fin (d N)) ℝ))
    [∀ k, IsProbabilityMeasure (L k)]
    (hlaw : HasLaw (fun ω i => (m.tbl i).Z N ω) (Measure.pi L) (μ N))
    (A : Matrix (Fin (n i N)) (Fin (d N)) ℝ → (Fin (d N) → ℝ)) (hAmeas : Measurable A)
    (hAunit : ∀ Z, A Z ⬝ᵥ A Z = 1)
    (F : (Fin (d N) → ℝ) → Matrix (Fin (n j N)) (Fin (d N)) ℝ → ℝ)
    (hFmeas : Measurable fun q : (Fin (d N) → ℝ) × Matrix (Fin (n j N)) (Fin (d N)) ℝ =>
      F q.1 q.2)
    (hFint : ∀ a : Fin (d N) → ℝ, a ⬝ᵥ a = 1 → Integrable (fun Z => F a Z ^ 2) (L j))
    {C ε : ℝ} (hε : 0 < ε)
    (hC : ∀ a : Fin (d N) → ℝ, a ⬝ᵥ a = 1 → ∫ Z, F a Z ^ 2 ∂(L j) ≤ C) :
    μ N {ω | ε ≤ |F (A ((m.tbl i).Z N ω)) ((m.tbl j).Z N ω)|}
      ≤ ENNReal.ofReal (C / ε ^ 2) := by
  classical
  have hprob : ∀ k : Fin M, IsProbabilityMeasure (L k) := fun k => inferInstance
  have hSmeas : MeasurableSet {q : Matrix (Fin (n j N)) (Fin (d N)) ℝ ×
      Matrix (Fin (n i N)) (Fin (d N)) ℝ | ε ≤ |F (A q.2) q.1|} := by
    have h1 : Measurable fun q : Matrix (Fin (n j N)) (Fin (d N)) ℝ ×
        Matrix (Fin (n i N)) (Fin (d N)) ℝ => F (A q.2) q.1 :=
      hFmeas.comp ((hAmeas.comp measurable_snd).prodMk measurable_fst)
    exact measurableSet_le measurable_const h1.abs
  have hpairmeas : MeasurableSet {Zs : (k : Fin M) → Matrix (Fin (n k N)) (Fin (d N)) ℝ |
      (Zs j, Zs i) ∈ {q : Matrix (Fin (n j N)) (Fin (d N)) ℝ ×
        Matrix (Fin (n i N)) (Fin (d N)) ℝ | ε ≤ |F (A q.2) q.1|}} :=
    ((measurable_pi_apply j).prodMk (measurable_pi_apply i)) hSmeas
  have key : μ N {ω | ε ≤ |F (A ((m.tbl i).Z N ω)) ((m.tbl j).Z N ω)|}
      = (Measure.pi L)
        {Zs | (Zs j, Zs i) ∈ {q : Matrix (Fin (n j N)) (Fin (d N)) ℝ ×
          Matrix (Fin (n i N)) (Fin (d N)) ℝ | ε ≤ |F (A q.2) q.1|}} :=
    hlaw.measure_eq hpairmeas
  rw [key]
  refine pi_pair_le L hprob hij.symm hSmeas _ ?_
  intro zi
  have hmeas' : Measurable fun zj : Matrix (Fin (n j N)) (Fin (d N)) ℝ => |F (A zi) zj| :=
    (hFmeas.comp (measurable_const.prodMk measurable_id)).abs
  have hint' : Integrable (fun zj => |F (A zi) zj| ^ 2) (L j) := by
    simpa [sq_abs] using hFint (A zi) (hAunit zi)
  have hC' : ∫ zj, |F (A zi) zj| ^ 2 ∂(L j) ≤ C := by
    simpa [sq_abs] using hC (A zi) (hAunit zi)
  exact R2.meas_ge_le_of_integral_sq _ hmeas' hint' hε hC'

/-- Item P with the direction read off table `i`, general noise: `‖E_j v̂_i‖² → c_j`. -/
theorem noise_projection_topDir_tendsto_general (m : MultiTableModel μ M n d) {ν : Measure ℝ}
    (hν : NoiseLaw ν) (hG : m.JointGeneralNoise ν) {i j : Fin M} (hij : i ≠ j) {cj : ℝ}
    (hregj : (m.tbl j).Regime cj) :
    TendstoInProb μ (fun N ω => ThetaEst.sqNormMulVec ((m.tbl j).E N ω)
      (ThetaEst.topDir ((m.tbl i).X N ω))) cj := by
  have := hν.prob
  have hd : Tendsto d atTop atTop := hregj.2.1
  have hbnd : ∀ ε : ℝ, 0 < ε → ∀ N,
      μ N {ω | ε ≤ |ThetaEst.sqNormMulVec ((m.tbl j).E N ω)
          (ThetaEst.topDir ((m.tbl i).X N ω)) - (n j N : ℝ) / d N|}
        ≤ ENNReal.ofReal (((∫ x, x ^ 4 ∂ν) + 2) * (n j N) / (d N : ℝ) ^ 2 / ε ^ 2) := by
    intro ε hε N
    have hD : 0 < d N := (m.tbl j).hd N
    set A : Matrix (Fin (n i N)) (Fin (d N)) ℝ → (Fin (d N) → ℝ) :=
      fun Z => WithLp.ofLp (ThetaEst.topDir (m.dataOf i N Z)) with hAdef
    set F : (Fin (d N) → ℝ) → Matrix (Fin (n j N)) (Fin (d N)) ℝ → ℝ :=
      fun a Z => ((d N : ℝ))⁻¹ * ((Z *ᵥ a) ⬝ᵥ (Z *ᵥ a)) - (n j N : ℝ) / d N with hFdef
    have hAmeas : Measurable A := by
      rw [hAdef]
      exact (WithLp.measurable_ofLp 2 (Fin (d N) → ℝ)).comp
        (ThetaEst.measurable_topDir.comp (m.measurable_dataOf i N))
    have hAunit : ∀ Z, A Z ⬝ᵥ A Z = 1 := by
      intro Z
      rw [hAdef]
      exact ThetaEst.dotProduct_topDir_self hD _
    have hFmeas : Measurable fun q : (Fin (d N) → ℝ) ×
        Matrix (Fin (n j N)) (Fin (d N)) ℝ => F q.1 q.2 := by
      rw [hFdef]
      exact (ThetaEst.measurable_dotProduct_pair.const_mul _).sub measurable_const
    have hFint : ∀ a : Fin (d N) → ℝ, a ⬝ᵥ a = 1 →
        Integrable (fun Z => F a Z ^ 2) (noiseMatrix ν (n j N) (d N)) := by
      intro a hau
      rw [hFdef]
      exact NoiseLaw.integrable_sq_normSqMulVec (r := n j N) hν hau
    have hset : {ω | ε ≤ |ThetaEst.sqNormMulVec ((m.tbl j).E N ω)
          (ThetaEst.topDir ((m.tbl i).X N ω)) - (n j N : ℝ) / d N|}
        = {ω | ε ≤ |F (A ((m.tbl i).Z N ω)) ((m.tbl j).Z N ω)|} := by
      ext ω
      rw [Set.mem_ofPred_eq, Set.mem_ofPred_eq, hFdef, hAdef]
      simp only
      rw [(m.tbl j).sqNormMulVec_E N ω (ThetaEst.topDir ((m.tbl i).X N ω)),
        m.dataOf_eq i N ω]
    rw [hset]
    refine m.measure_randomDir_le_general hij N (fun k => noiseMatrix ν (n k N) (d N))
      (hG N) A hAmeas hAunit F hFmeas hFint hε ?_
    intro a hau
    rw [hFdef]
    exact NoiseLaw.integral_sq_normSqMulVec_le (r := n j N) hν hau
  have hK : Tendsto (fun N => ((∫ x, x ^ 4 ∂ν) + 2) * (n j N : ℝ) / (d N : ℝ) ^ 2)
      atTop (𝓝 0) := by
    have hdR : Tendsto (fun N => ((d N : ℝ))) atTop atTop :=
      tendsto_natCast_atTop_atTop.comp hd
    have hinv : Tendsto (fun N => ((d N : ℝ))⁻¹) atTop (𝓝 0) := hdR.inv_tendsto_atTop
    have hprod := (hregj.2.2.const_mul ((∫ x, x ^ 4 ∂ν) + 2)).mul hinv
    have hz : ((∫ x, x ^ 4 ∂ν) + 2) * cj * 0 = 0 := by ring
    rw [hz] at hprod
    exact hprod.congr fun N => by ring
  have h0 := ThetaEst.tendstoInProb_of_meas_le hK hbnd
  have hfin := h0.add_tendsto hregj.2.2
  simpa using hfin

/-- The cross term `u_jᵀ E_j v̂_i → 0`, general noise. -/
theorem cross_term_tendsto_general (m : MultiTableModel μ M n d) {ν : Measure ℝ}
    (hν : NoiseLaw ν) (hG : m.JointGeneralNoise ν) {i j : Fin M} (hij : i ≠ j)
    (hd : Tendsto d atTop atTop) :
    TendstoInProb μ (fun N ω => WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ
      ((m.tbl j).E N ω *ᵥ WithLp.ofLp (ThetaEst.topDir ((m.tbl i).X N ω)))) 0 := by
  have := hν.prob
  have hbnd : ∀ ε : ℝ, 0 < ε → ∀ N,
      μ N {ω | ε ≤ |WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ
          ((m.tbl j).E N ω *ᵥ WithLp.ofLp (ThetaEst.topDir ((m.tbl i).X N ω)))
          - (fun _ : ℕ => (0 : ℝ)) N|}
        ≤ ENNReal.ofReal (((d N : ℝ))⁻¹ / ε ^ 2) := by
    intro ε hε N
    have hD : 0 < d N := (m.tbl j).hd N
    have hxu : WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ WithLp.ofLp ((m.tbl j).u N) = 1 :=
      (m.tbl j).dotProduct_u_self N
    have hs2 : ((Real.sqrt (d N))⁻¹) ^ 2 = ((d N : ℝ))⁻¹ := by
      rw [sq]
      exact R2.sqrt_inv_mul_sqrt_inv (d N)
    set A : Matrix (Fin (n i N)) (Fin (d N)) ℝ → (Fin (d N) → ℝ) :=
      fun Z => WithLp.ofLp (ThetaEst.topDir (m.dataOf i N Z)) with hAdef
    set F : (Fin (d N) → ℝ) → Matrix (Fin (n j N)) (Fin (d N)) ℝ → ℝ :=
      fun a Z => (Real.sqrt (d N))⁻¹ * (WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ (Z *ᵥ a)) with hFdef
    have hsq : ∀ (a : Fin (d N) → ℝ) (Z : Matrix (Fin (n j N)) (Fin (d N)) ℝ),
        F a Z ^ 2
          = ((d N : ℝ))⁻¹ * (WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ (Z *ᵥ a)) ^ 2 := by
      intro a Z
      rw [hFdef]
      simp only
      rw [mul_pow, hs2]
    have hAmeas : Measurable A := by
      rw [hAdef]
      exact (WithLp.measurable_ofLp 2 (Fin (d N) → ℝ)).comp
        (ThetaEst.measurable_topDir.comp (m.measurable_dataOf i N))
    have hAunit : ∀ Z, A Z ⬝ᵥ A Z = 1 := by
      intro Z
      rw [hAdef]
      exact ThetaEst.dotProduct_topDir_self hD _
    have hFmeas : Measurable fun q : (Fin (d N) → ℝ) ×
        Matrix (Fin (n j N)) (Fin (d N)) ℝ => F q.1 q.2 := by
      rw [hFdef]
      exact (ThetaEst.measurable_dotProduct_left_pair
        (WithLp.ofLp ((m.tbl j).u N))).const_mul _
    have hFint : ∀ a : Fin (d N) → ℝ, a ⬝ᵥ a = 1 →
        Integrable (fun Z => F a Z ^ 2) (noiseMatrix ν (n j N) (d N)) := by
      intro a _
      simp only [hsq]
      exact (NoiseLaw.integrable_sq_dotProduct_mulVec (r := n j N) hν
        (WithLp.ofLp ((m.tbl j).u N)) a).const_mul _
    have hset : {ω | ε ≤ |WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ
          ((m.tbl j).E N ω *ᵥ WithLp.ofLp (ThetaEst.topDir ((m.tbl i).X N ω)))
          - (fun _ : ℕ => (0 : ℝ)) N|}
        = {ω | ε ≤ |F (A ((m.tbl i).Z N ω)) ((m.tbl j).Z N ω)|} := by
      ext ω
      rw [Set.mem_ofPred_eq, Set.mem_ofPred_eq, sub_zero, hFdef, hAdef]
      simp only
      have hE : (m.tbl j).E N ω = (Real.sqrt (d N))⁻¹ • (m.tbl j).Z N ω := rfl
      rw [hE, Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul, m.dataOf_eq i N ω]
    rw [hset]
    refine m.measure_randomDir_le_general hij N (fun k => noiseMatrix ν (n k N) (d N))
      (hG N) A hAmeas hAunit F hFmeas hFint hε ?_
    intro a hau
    simp only [hsq]
    rw [integral_const_mul,
      NoiseLaw.integral_sq_dotProduct_mulVec (r := n j N) hν hxu hau, mul_one]
  have hK : Tendsto (fun N => ((d N : ℝ))⁻¹) atTop (𝓝 0) := by
    have hdR : Tendsto (fun N => ((d N : ℝ))) atTop atTop :=
      tendsto_natCast_atTop_atTop.comp hd
    exact hdR.inv_tendsto_atTop
  have h0 := ThetaEst.tendstoInProb_of_meas_le (ctr := fun _ => (0 : ℝ)) hK hbnd
  simpa using h0

/-- The general-noise discharge of `ThetaEstLaw`: independent tables `i ≠ j` with i.i.d. entries
of a fixed law `ν` of mean `0`, variance `1` and finite fourth moment, in the regime of table
`j`, give both limits. -/
theorem thetaEstLaw_of_general (m : MultiTableModel μ M n d) {ν : Measure ℝ} (hν : NoiseLaw ν)
    {i j : Fin M} (hij : i ≠ j) (hG : m.JointGeneralNoise ν) {cj : ℝ}
    (hregj : (m.tbl j).Regime cj) : m.ThetaEstLaw i j cj :=
  ⟨m.noise_projection_topDir_tendsto_general hν hG hij hregj,
    m.cross_term_tendsto_general hν hG hij hregj.2.1⟩

/-- The Gaussian case is the instance `ν = gaussianReal 0 1` (`noiseLaw_gaussian`); this checks
that the general theorem gives `thetaEstLaw_of_gaussian` back. -/
example (m : MultiTableModel μ M n d) {i j : Fin M} (hij : i ≠ j) (hG : m.JointGaussianNoise)
    {cj : ℝ} (hregj : (m.tbl j).Regime cj) : m.ThetaEstLaw i j cj :=
  m.thetaEstLaw_of_general noiseLaw_gaussian hij hG hregj

/-- **`thm:theta_est` at a general noise law** (`assum:general_noise`, `main_paper.tex:244`):
the estimator is consistent for `θ_j` when the noise entries are i.i.d. with mean `0`,
variance `1` and a finite fourth moment. `SingleTableLaw` of table `i` stays a hypothesis: its
general-noise discharge is a later stage. -/
theorem thm_theta_est_general [∀ N, IsProbabilityMeasure (μ N)] (m : MultiTableModel μ M n d)
    {ν : Measure ℝ} (hν : NoiseLaw ν) {i j : Fin M} (hij : i ≠ j) {ci cj : ℝ} (hci : 0 ≤ ci)
    (hthr : ci < (m.tbl i).θ ^ 4) (law : (m.tbl i).SingleTableLaw ci)
    (hG : m.JointGeneralNoise ν) (hregj : (m.tbl j).Regime cj) :
    TendstoInProb μ (fun N ω => m.thetaHat i j ci cj N ω) (m.tbl j).θ :=
  m.thm_theta_est hci hthr law (m.thetaEstLaw_of_general hν hij hG hregj)

end MultiTableModel

end StackedSVD

