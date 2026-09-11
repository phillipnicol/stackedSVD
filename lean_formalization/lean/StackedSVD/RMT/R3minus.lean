/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.R5
import StackedSVD.RMT.R1
import StackedSVD.RMT.T

/-!
# Item R3⁻: the lower edge of the Wishart block

Review note: `notes/archive/rmt_R3.md`, section "Statement (lower edge, item R3⁻)". Item R3 says the
top eigenvalue of `W₀` does not exceed `bulkEdge c`; this file says it does not fall below it:
for every `δ > 0`,

`Tendsto (fun N => μ N {ω | lamMax (W₀ N ω) _ ≤ bulkEdge c - δ}) atTop (𝓝 0)`.

The argument has four steps, in the order of the note.

1. `trace2_le_of_lamMax_le`: deterministic. If `lamMax W ≤ b - δ` and `b ≤ x` then every
   eigenvalue satisfies `x - λ_a ≥ δ`, so `d⁻¹ tr G(x)² ≤ (δ²)⁻¹`.
2. `exists_mDeriv_gt`: from `MP.mDeriv_tendsto_atTop`, every level `K` is exceeded by
   `MP.mDeriv c x₀` at some real `x₀ ∈ (b, b + 1)`.
3. `tendstoInProb_trace2_of_complex`: at a real `x > b` the trace `d⁻¹ tr G(x)²` converges in
   probability to `MP.mDeriv c x`. Inputs: R1b at `x + iη` (`R1.tendstoInProb_stieltjes2C`),
   the transfer bound T2t on the edge event (`T.norm_trace2_sub_stieltjes2C_le`, error
   `16η/ε³`) and `MP.tendsto_mCDeriv` for `η ↓ 0`. The `ε`-`η` bookkeeping mirrors
   `T.tendstoInProb_cform2_of_complex'`, with two bad sets instead of three: the trace needs no
   bound on a test vector.
4. `tendsto_measure_lamMax_le_sub`: on the lower-edge event step 1 caps the trace at `(δ²)⁻¹`,
   while step 3 puts it near `MP.mDeriv c x₀ > 2 (δ²)⁻¹` (step 2 at `K := 2 (δ²)⁻¹`), so the
   event sits inside `{(δ²)⁻¹ ≤ |trace - MP.mDeriv c x₀|}`, whose measure tends to `0`.

The `SpikedModel` section then instantiates all four at the model of `RMT/R0.lean`: the
complex input comes from R0's Gaussian block and R1b, the Gram version follows from
`R4.lamMax_le_lamMax_vecMulVec` with `R0.gram_eq`, and, with R5's
`tendsto_measure_lamMax_le_of_subcritical`, the two bounds give the `lamMax` field of
`SingleTableLaw` in the subcritical regime `θ⁴ ≤ c`, where `rhoSq θ c = bulkEdge c`.

STATUS 2026-08-30: `lean -j 3 -R . StackedSVD/RMT/R3minus.lean` exit 0, 0 `sorry`, no warning.
The five `R3minus` declarations use only the standard axioms; the five `SpikedModel` ones
inherit `sorryAx` from the three open `sorry` of `RMT/R1.lean` (closed 2026-08-30: R1 is
complete, see `RMT/SteinStep.lean`), which another task closes.
See `notes/archive/agent_reports/proof_r3minus.md`.
-/

open Filter Topology MeasureTheory ProbabilityTheory
open scoped Matrix ENNReal

namespace StackedSVD
namespace R3minus

/-! ### Step 1: the deterministic cap on the trace -/

section Deterministic

variable {d : ℕ} {W : Matrix (Fin d) (Fin d) ℝ}

/-- The trace of `G₀(x)²` as an eigenvalue sum, above the top eigenvalue. This is the public
copy of the private `T.trace_resolv_sq_eq_sum`. -/
theorem trace_resolv_sq_eq_sum (hW : W.IsHermitian) {x : ℝ} (hlt : lamMax W hW < x) :
    (R4.resolv W x * R4.resolv W x).trace = ∑ a, ((hW.eigenvalues a - x)⁻¹) ^ 2 := by
  rw [R4.resolv_mul_resolv_eq_conj hW hlt, Matrix.trace_mul_comm, ← Matrix.mul_assoc,
    R4.transpose_eigU_mul, Matrix.one_mul, Matrix.trace_diagonal]

/-- **Step 1.** On the lower-edge event `lamMax W ≤ b - δ` the trace of the squared resolvent
stays below `(δ²)⁻¹` at every real `x` at or above `b`: every eigenvalue is at distance at
least `δ` from `x`. -/
theorem trace2_le_of_lamMax_le (hW : W.IsHermitian) {b δ x : ℝ} (hδ : 0 < δ)
    (hlam : lamMax W hW ≤ b - δ) (hx : b ≤ x) :
    (d : ℝ)⁻¹ * (R4.resolv W x * R4.resolv W x).trace ≤ (δ ^ 2)⁻¹ := by
  have hlt : lamMax W hW < x := by linarith
  rw [trace_resolv_sq_eq_sum hW hlt]
  have hterm : ∀ a : Fin d, ((hW.eigenvalues a - x)⁻¹) ^ 2 ≤ (δ ^ 2)⁻¹ := by
    intro a
    have h1 : hW.eigenvalues a ≤ b - δ := le_trans (R4.eigenvalues_le_lamMax hW a) hlam
    have h2 : hW.eigenvalues a - x ≤ -δ := by linarith
    have h3 : δ ^ 2 ≤ (hW.eigenvalues a - x) ^ 2 := by nlinarith
    calc ((hW.eigenvalues a - x)⁻¹) ^ 2 = ((hW.eigenvalues a - x) ^ 2)⁻¹ := by rw [inv_pow]
      _ ≤ (δ ^ 2)⁻¹ := inv_anti₀ (by positivity) h3
  have hsum : ∑ a : Fin d, ((hW.eigenvalues a - x)⁻¹) ^ 2 ≤ (d : ℝ) * (δ ^ 2)⁻¹ := by
    calc ∑ a : Fin d, ((hW.eigenvalues a - x)⁻¹) ^ 2 ≤ ∑ _a : Fin d, (δ ^ 2)⁻¹ :=
          Finset.sum_le_sum fun a _ => hterm a
      _ = (d : ℝ) * (δ ^ 2)⁻¹ := by
          rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  rcases Nat.eq_zero_or_pos d with hd | hd
  · subst hd
    simp only [Finset.univ_eq_empty, Finset.sum_empty, mul_zero]
    positivity
  · have hdR : (0 : ℝ) < d := by exact_mod_cast hd
    calc (d : ℝ)⁻¹ * ∑ a : Fin d, ((hW.eigenvalues a - x)⁻¹) ^ 2
        ≤ (d : ℝ)⁻¹ * ((d : ℝ) * (δ ^ 2)⁻¹) := by
          exact mul_le_mul_of_nonneg_left hsum (by positivity)
      _ = (δ ^ 2)⁻¹ := by field_simp

end Deterministic

/-! ### Step 2: `MP.mDeriv` exceeds every level just above the edge -/

/-- **Step 2.** `MP.mDeriv c` tends to `+∞` as `x` decreases to `bulkEdge c`, so every level
`K` is exceeded at some real point of `(bulkEdge c, bulkEdge c + 1)`. -/
theorem exists_mDeriv_gt {c : ℝ} (hc : 0 < c) (K : ℝ) :
    ∃ x : ℝ, bulkEdge c < x ∧ x < bulkEdge c + 1 ∧ K < MP.mDeriv c x := by
  have hK : ∀ᶠ x : ℝ in 𝓝[>] (bulkEdge c), K < MP.mDeriv c x :=
    (MP.mDeriv_tendsto_atTop hc).eventually_gt_atTop K
  have hsmall : ∀ᶠ x : ℝ in 𝓝[>] (bulkEdge c), x < bulkEdge c + 1 :=
    Filter.Eventually.filter_mono nhdsWithin_le_nhds
      (Filter.eventually_of_mem (Iio_mem_nhds (by linarith)) fun _ hy => hy)
  have hev : ∀ᶠ x : ℝ in 𝓝[>] (bulkEdge c),
      bulkEdge c < x ∧ x < bulkEdge c + 1 ∧ K < MP.mDeriv c x := by
    filter_upwards [self_mem_nhdsWithin, hsmall, hK] with x h1 h2 h3
    exact ⟨h1, h2, h3⟩
  exact hev.exists

/-! ### Steps 3 and 4: the probabilistic half -/

section Prob

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {dN : ℕ → ℕ}
  {c : ℝ}

/-- **Step 3.** At a real `x > bulkEdge c` the trace `d⁻¹ tr G₀(x)²` converges in probability
to `MP.mDeriv c x`. `hcplx` is R1b in its norm form, at every complex point of the upper half
plane; `hedge` is item R3 in complement form. -/
theorem tendstoInProb_trace2_of_complex (hc : 0 < c)
    (W₀ : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ)
    (hsymm : ∀ N ω, (W₀ N ω).IsHermitian) (hdpos : ∀ N, 0 < dN N)
    (hedge : ∀ ε > 0, Tendsto
      (fun N => μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ bulkEdge c + ε}ᶜ) atTop (𝓝 0))
    (hcplx : ∀ z : ℂ, 0 < z.im →
      TendstoInProb μ (fun N ω => ‖R4C.stieltjes2C (W₀ N ω) z - MP.mCDeriv c z‖) 0)
    {x : ℝ} (hx : bulkEdge c < x) :
    TendstoInProb μ (fun N ω =>
      (dN N : ℝ)⁻¹ * (R4.resolv (W₀ N ω) x * R4.resolv (W₀ N ω) x).trace) (MP.mDeriv c x) := by
  refine tendstoInProb_of_subset_union₂ fun δ hδ => ?_
  set ε : ℝ := x - bulkEdge c with hεdef
  have hε : 0 < ε := by rw [hεdef]; linarith
  have hxε : bulkEdge c + ε ≤ x := by rw [hεdef]; linarith
  -- the choice of `η`: small enough for the T2t error and for the scalar limit
  obtain ⟨η, hη, hηsmall, hηL⟩ : ∃ η : ℝ, 0 < η ∧ η < δ * ε ^ 3 / 192 ∧
      ‖MP.mCDeriv c ((x : ℂ) + (η : ℂ) * Complex.I) - ((MP.mDeriv c x : ℝ) : ℂ)‖ < δ / 3 := by
    have hball : Metric.ball (((MP.mDeriv c x : ℝ)) : ℂ) (δ / 3)
        ∈ 𝓝 (((MP.mDeriv c x : ℝ)) : ℂ) := Metric.ball_mem_nhds _ (by positivity)
    have hLnear : ∀ᶠ η : ℝ in 𝓝[>] (0 : ℝ),
        ‖MP.mCDeriv c ((x : ℂ) + (η : ℂ) * Complex.I) - ((MP.mDeriv c x : ℝ) : ℂ)‖ < δ / 3 := by
      filter_upwards [MP.tendsto_mCDeriv hc hx hball] with η hη
      simpa [Metric.mem_ball, dist_eq_norm] using hη
    have hr : (0 : ℝ) < δ * ε ^ 3 / 192 := by positivity
    have hsmall : ∀ᶠ η : ℝ in 𝓝[>] (0 : ℝ), η < δ * ε ^ 3 / 192 :=
      Filter.Eventually.filter_mono nhdsWithin_le_nhds
        (Filter.eventually_of_mem (Iio_mem_nhds hr) fun _ hy => hy)
    have hev : ∀ᶠ η : ℝ in 𝓝[>] (0 : ℝ), 0 < η ∧ η < δ * ε ^ 3 / 192 ∧
        ‖MP.mCDeriv c ((x : ℂ) + (η : ℂ) * Complex.I) - ((MP.mDeriv c x : ℝ) : ℂ)‖ < δ / 3 := by
      filter_upwards [self_mem_nhdsWithin, hsmall, hLnear] with η h1 h2 h3
      exact ⟨h1, h2, h3⟩
    exact hev.exists
  have hzim : (0 : ℝ) < ((x : ℂ) + (η : ℂ) * Complex.I).im := by simpa using hη
  refine ⟨fun N => {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ bulkEdge c + ε / 2}ᶜ,
    fun N => {ω | δ / 3 ≤ |‖R4C.stieltjes2C (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I)
      - MP.mCDeriv c ((x : ℂ) + (η : ℂ) * Complex.I)‖ - 0|},
    ?_, hedge (ε / 2) (by positivity),
    hcplx ((x : ℂ) + (η : ℂ) * Complex.I) hzim (δ / 3) (by positivity)⟩
  intro N ω hω
  by_contra hnot
  have hlam : lamMax (W₀ N ω) (hsymm N ω) ≤ bulkEdge c + ε / 2 := by
    by_contra hcon
    exact hnot (Set.mem_union_left _ (Set.mem_compl hcon))
  have hCf : ‖R4C.stieltjes2C (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I)
      - MP.mCDeriv c ((x : ℂ) + (η : ℂ) * Complex.I)‖ < δ / 3 := by
    by_contra hcon
    refine hnot (Set.mem_union_right _ ?_)
    change δ / 3 ≤ |‖R4C.stieltjes2C (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I)
      - MP.mCDeriv c ((x : ℂ) + (η : ℂ) * Complex.I)‖ - 0|
    rw [sub_zero, abs_norm]
    exact not_lt.mp hcon
  have hT2 := T.norm_trace2_sub_stieltjes2C_le (W := W₀ N ω) (hsymm N ω) hc hε hη
    (hdpos N) hlam hxε
  have hstep : 16 * η / ε ^ 3 < δ / 3 := by
    rw [div_lt_iff₀ (by positivity)]
    nlinarith [hηsmall, pow_pos hε 3]
  -- the triangle inequality on the three gaps
  set Tr : ℝ := (dN N : ℝ)⁻¹ * (R4.resolv (W₀ N ω) x * R4.resolv (W₀ N ω) x).trace with hTrdef
  have habs : |Tr - MP.mDeriv c x| = ‖((Tr : ℝ) : ℂ) - ((MP.mDeriv c x : ℝ) : ℂ)‖ := by
    rw [← Complex.ofReal_sub, Complex.norm_real, Real.norm_eq_abs]
  have hδle : δ ≤ |Tr - MP.mDeriv c x| := hω
  rw [habs] at hδle
  have htri : ‖((Tr : ℝ) : ℂ) - ((MP.mDeriv c x : ℝ) : ℂ)‖
      ≤ ‖((Tr : ℝ) : ℂ) - R4C.stieltjes2C (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I)‖
        + (‖R4C.stieltjes2C (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I)
            - MP.mCDeriv c ((x : ℂ) + (η : ℂ) * Complex.I)‖
          + ‖MP.mCDeriv c ((x : ℂ) + (η : ℂ) * Complex.I) - ((MP.mDeriv c x : ℝ) : ℂ)‖) := by
    have e : ((Tr : ℝ) : ℂ) - ((MP.mDeriv c x : ℝ) : ℂ)
        = (((Tr : ℝ) : ℂ) - R4C.stieltjes2C (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I))
          + ((R4C.stieltjes2C (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I)
              - MP.mCDeriv c ((x : ℂ) + (η : ℂ) * Complex.I))
            + (MP.mCDeriv c ((x : ℂ) + (η : ℂ) * Complex.I)
              - ((MP.mDeriv c x : ℝ) : ℂ))) := by ring
    rw [e]
    exact (norm_add_le _ _).trans (add_le_add le_rfl (norm_add_le _ _))
  linarith [hT2, hstep, hCf, hηL, hδle, htri]

/-- **Step 4, R3⁻.** The lower edge. `htr` is the trace form of (H2) at real points, which
step 3 builds from R1b and item T. -/
theorem tendsto_measure_lamMax_le_sub (hc : 0 < c)
    (W₀ : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ)
    (hsymm : ∀ N ω, (W₀ N ω).IsHermitian)
    (htr : ∀ x : ℝ, bulkEdge c < x → TendstoInProb μ (fun N ω =>
      (dN N : ℝ)⁻¹ * (R4.resolv (W₀ N ω) x * R4.resolv (W₀ N ω) x).trace) (MP.mDeriv c x))
    {δ : ℝ} (hδ : 0 < δ) :
    Tendsto (fun N => μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ bulkEdge c - δ}) atTop (𝓝 0) := by
  have hpos : (0 : ℝ) < (δ ^ 2)⁻¹ := by positivity
  obtain ⟨x₀, hx₀lo, -, hx₀K⟩ := exists_mDeriv_gt hc (2 * (δ ^ 2)⁻¹)
  have hlim := htr x₀ hx₀lo (δ ^ 2)⁻¹ hpos
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hlim
    (fun _ => zero_le) fun N => measure_mono ?_
  intro ω hω
  have hcap : (dN N : ℝ)⁻¹ * (R4.resolv (W₀ N ω) x₀ * R4.resolv (W₀ N ω) x₀).trace
      ≤ (δ ^ 2)⁻¹ := trace2_le_of_lamMax_le (hsymm N ω) hδ hω hx₀lo.le
  change (δ ^ 2)⁻¹ ≤ |(dN N : ℝ)⁻¹ * (R4.resolv (W₀ N ω) x₀ * R4.resolv (W₀ N ω) x₀).trace
    - MP.mDeriv c x₀|
  rw [abs_sub_comm, abs_of_nonneg (by linarith)]
  linarith

end Prob

end R3minus

/-! ### The model: `W₀`, the Gram matrix, and the subcritical `lamMax` field -/

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- `ω ↦ X N ω` is measurable, entrywise, from `m.hZ`. Same route as `measurable_Eperp` of
`RMT/R5.lean`. -/
theorem measurable_X (m : SpikedModel μ n d) (N : ℕ) : Measurable (m.X N) := by
  have hZ : ∀ (r : Fin (n N)) (j : Fin (d N)), Measurable fun ω => m.Z N ω r j := by
    intro r j
    have h1 : Measurable fun ω => m.Z N ω r := (measurable_pi_apply r).comp (m.hZ N)
    exact (measurable_pi_apply j).comp h1
  refine measurable_pi_lambda _ fun r => measurable_pi_lambda _ fun j => ?_
  have h : (fun ω => m.X N ω r j)
      = fun ω => m.θ * (WithLp.ofLp (m.u N) r * WithLp.ofLp (m.v N) j)
        + (Real.sqrt (d N))⁻¹ * m.Z N ω r j := rfl
  rw [h]
  exact measurable_const.add ((hZ r j).const_mul _)

/-- The upper level sets of `gramLamMax` are measurable. R5 proves the companion
`measurableSet_lamMax_W0_le` for `W₀`. -/
theorem measurableSet_gramLamMax_le (m : SpikedModel μ n d) (N : ℕ) (r : ℝ) :
    MeasurableSet {ω | gramLamMax (m.X N ω) ≤ r} := by
  have h : {ω | gramLamMax (m.X N ω) ≤ r}
      = (fun ω => gramLamMax (m.X N ω)) ⁻¹' Set.Iic r := rfl
  rw [h]
  exact (measurable_gramLamMax.comp (m.measurable_X N)) measurableSet_Iic

/-- The rank-one part only raises the top eigenvalue: `λ_max(W₀) ≤ λ_max(Xᵀ X)`. This is
`R4.lamMax_le_lamMax_vecMulVec` transported along R0's `gram_eq`. -/
theorem lamMax_W0_le_gramLamMax (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) :
    lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ gramLamMax (m.X N ω) := by
  have hgram := m.gram_eq N ω
  have hH : (m.W0 N ω + Matrix.vecMulVec (m.qvec N ω) (m.qvec N ω)).IsHermitian := by
    rw [← hgram]; exact isHermitian_transpose_mul_self _
  calc lamMax (m.W0 N ω) (m.isHermitian_W0 N ω)
      ≤ lamMax (m.W0 N ω + Matrix.vecMulVec (m.qvec N ω) (m.qvec N ω)) hH :=
        R4.lamMax_le_lamMax_vecMulVec (m.isHermitian_W0 N ω) hH
    _ = gramLamMax (m.X N ω) :=
        (lamMax_congr hgram (isHermitian_transpose_mul_self (m.X N ω)) hH).symm

variable [∀ N, IsProbabilityMeasure (μ N)] {m : SpikedModel μ n d} {c : ℝ}

omit [∀ N, IsProbabilityMeasure (μ N)] in
/-- R1b on the model's `W₀`. R0's block `B ~ gaussianMatrix (n - 1) d` with
`W₀ = d⁻¹ Bᵀ B` (`exists_block_hasLaw_snd`) turns `R1.tendstoInProb_stieltjes2C` into a
statement about `m.W0`. The hypothesis `2 ≤ n N` is what makes the block nonempty; see the
statement choices of `notes/archive/agent_reports/proof_r3minus.md`. -/
theorem tendstoInProb_stieltjes2C_W0 (hc : 0 < c) (hreg : m.Regime c) (hG : m.GaussianNoise)
    (hn2 : ∀ N, 2 ≤ n N) {z : ℂ} (hz : 0 < z.im) :
    TendstoInProb μ (fun N ω => ‖R4C.stieltjes2C (m.W0 N ω) z - MP.mCDeriv c z‖) 0 := by
  classical
  have hpsucc : ∀ N, n N = (n N - 1) + 1 := fun N =>
    (Nat.succ_pred_eq_of_pos (m.hn N)).symm
  choose B hW hlaw using fun N => m.exists_block_hasLaw_snd N hG (hpsucc N)
  have hppos : ∀ N, 0 < n N - 1 := fun N => by have := hn2 N; omega
  have hR1 := R1.tendstoInProb_stieltjes2C (c := c) (z := z) (pN := fun N => n N - 1)
    (dN := d) hc hz hreg.2.1 hppos (m.tendsto_pred_ratio hreg) B hlaw
  have heq : (fun N (ω : Ω N) => ‖R4C.stieltjes2C (R1.W0 (B N ω)) z - MP.mCDeriv c z‖)
      = fun N (ω : Ω N) => ‖R4C.stieltjes2C (m.W0 N ω) z - MP.mCDeriv c z‖ := by
    funext N ω
    have : R1.W0 (B N ω) = m.W0 N ω := (hW N ω).symm
    rw [this]
  rwa [heq] at hR1

/-- The trace form of (H2) at a real `x > bulkEdge c`, for the model's `W₀`: step 3 at the
edge bound of `ResolventLimits` and at R1b through `tendstoInProb_stieltjes2C_W0`. -/
theorem tendstoInProb_trace2_W0 (H : m.ResolventLimits c) (hc : 0 < c) (hreg : m.Regime c)
    (hG : m.GaussianNoise) (hn2 : ∀ N, 2 ≤ n N) {x : ℝ} (hx : bulkEdge c < x) :
    TendstoInProb μ (fun N ω =>
      (d N : ℝ)⁻¹ * (R4.resolv (m.W0 N ω) x * R4.resolv (m.W0 N ω) x).trace)
      (MP.mDeriv c x) := by
  refine R3minus.tendstoInProb_trace2_of_complex hc (fun N ω => m.W0 N ω)
    (fun N ω => m.isHermitian_W0 N ω) (fun N => m.hd N) (fun ε hε => ?_)
    (fun z hz => tendstoInProb_stieltjes2C_W0 hc hreg hG hn2 hz) hx
  exact tendsto_measure_compl_zero
    (fun N => (m.measurableSet_lamMax_W0_le N _).nullMeasurableSet) (H.edge ε hε)

/-- **R3⁻ for the model.** The top eigenvalue of `W₀` does not fall below `bulkEdge c`. -/
theorem tendsto_measure_lamMax_W0_le_sub (H : m.ResolventLimits c) (hc : 0 < c)
    (hreg : m.Regime c) (hG : m.GaussianNoise) (hn2 : ∀ N, 2 ≤ n N) {δ : ℝ} (hδ : 0 < δ) :
    Tendsto (fun N => μ N {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ bulkEdge c - δ})
      atTop (𝓝 0) :=
  R3minus.tendsto_measure_lamMax_le_sub hc (fun N ω => m.W0 N ω)
    (fun N ω => m.isHermitian_W0 N ω)
    (fun _ hx => tendstoInProb_trace2_W0 H hc hreg hG hn2 hx) hδ

/-- **R3⁻ for the Gram matrix.** `λ_max(Xᵀ X) ≥ λ_max(W₀)`, so the Gram top eigenvalue does
not fall below `bulkEdge c` either. -/
theorem tendsto_measure_gramLamMax_le_sub (H : m.ResolventLimits c) (hc : 0 < c)
    (hreg : m.Regime c) (hG : m.GaussianNoise) (hn2 : ∀ N, 2 ≤ n N) {δ : ℝ} (hδ : 0 < δ) :
    Tendsto (fun N => μ N {ω | gramLamMax (m.X N ω) ≤ bulkEdge c - δ}) atTop (𝓝 0) :=
  tendsto_measure_zero_of_subset
    (fun N _ hω => le_trans (m.lamMax_W0_le_gramLamMax N _) hω)
    (tendsto_measure_lamMax_W0_le_sub H hc hreg hG hn2 hδ)

/-- **The subcritical `lamMax` field of `SingleTableLaw`.** For `θ⁴ ≤ c` (the equality case
included) `rhoSq θ c = bulkEdge c`, and the two edges pin the Gram top eigenvalue there: R5's
`tendsto_measure_lamMax_le_of_subcritical` gives the upper bound, R3⁻ the lower one. -/
theorem lamMax_tendstoInProb_of_subcritical (H : m.ResolventLimits c) (hc : 0 < c)
    (hreg : m.Regime c) (hG : m.GaussianNoise) (hn2 : ∀ N, 2 ≤ n N) (hθ : m.θ ^ 4 ≤ c) :
    TendstoInProb μ (fun N ω => gramLamMax (m.X N ω)) (rhoSq m.θ c) := by
  have hrho : rhoSq m.θ c = bulkEdge c := by
    rw [rhoSq, if_neg (not_lt.mpr hθ)]
  rw [hrho]
  intro ε hε
  have hup : Tendsto (fun N => μ N {ω | bulkEdge c + ε ≤ gramLamMax (m.X N ω)}) atTop (𝓝 0) := by
    refine tendsto_measure_zero_of_subset (t := fun N =>
      {ω | gramLamMax (m.X N ω) ≤ bulkEdge c + ε / 2}ᶜ) (fun N ω hω => ?_) ?_
    · have hω' : bulkEdge c + ε ≤ gramLamMax (m.X N ω) := hω
      refine Set.mem_compl ?_
      intro hmem
      have h : gramLamMax (m.X N ω) ≤ bulkEdge c + ε / 2 := hmem
      linarith
    · exact tendsto_measure_compl_zero
        (fun N => (m.measurableSet_gramLamMax_le N _).nullMeasurableSet)
        (tendsto_measure_lamMax_le_of_subcritical H hc hθ (ε / 2) (by linarith))
  have hlo := tendsto_measure_gramLamMax_le_sub H hc hreg hG hn2 hε
  refine tendsto_measure_zero_of_subset (fun N ω hω => ?_)
    (tendsto_measure_zero_union hup hlo)
  have hω' : ε ≤ |gramLamMax (m.X N ω) - bulkEdge c| := hω
  rcases le_abs.mp hω' with h | h
  · exact Set.mem_union_left _ (show bulkEdge c + ε ≤ gramLamMax (m.X N ω) by linarith)
  · exact Set.mem_union_right _ (show gramLamMax (m.X N ω) ≤ bulkEdge c - ε by linarith)

end SpikedModel

end StackedSVD
