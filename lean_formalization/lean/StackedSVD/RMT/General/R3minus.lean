/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Trace
import StackedSVD.RMT.General.Companion
import StackedSVD.RMT.General.Stability
import StackedSVD.RMT.General.Defs
import StackedSVD.RMT.R5
import StackedSVD.RMT.T
import StackedSVD.RMT.MP7
import StackedSVD.Prob.TendstoInProb

/-!
# Item R3⁻ at a general noise law: the lower edge of the Wishart block

`notes/archive/prop_single_table_general.md` section 5, unit G9. General-law twin of
`RMT/R3minus.lean`. That file's `SpikedModel` section is Gaussian only because it routes R1b
through R0's block law `exists_block_hasLaw_snd`, which needs `2 ≤ n N` to make the block
nonempty. This file feeds unit G3's `tendstoInProb_stieltjes2C_general`
(`RMT/General/Trace.lean:1133`) instead, on `m.W0` through the downdate identity of unit G2
(`SpikedModel.W0_eq_gram_sub`, `RMT/General/Companion.lean:391`), so the five model-level
theorems copy across **without** `hn2`.

## Content

1. Namespace `GenRMT.R3minus`: the five model-free declarations of `RMT/R3minus.lean:62` to
   `:225`, copied verbatim except for the namespace. `RMT/General/` may not import
   `RMT/R3minus.lean` (choice 8 of the plan note: its closure holds `RMT/R1.lean`,
   `RMT/SteinStep.lean` and a `Vendor/COLT83/` file), so the four deterministic steps of item
   R3⁻ are duplicated here under `GenRMT.R3minus`, exactly as unit G3a duplicated the twelve
   stability lemmas of `RMT/R1.lean` into `RMT/General/Stability.lean`. Follow-up item F34 of
   `notes/FOLLOWUP_LIST.md` tracks the dedup.
2. Namespace `GenRMT`: `norm_stieltjesC_sub_vecMulVec_le` and `norm_stieltjes2C_sub_vecMulVec_le`,
   the general-`W`, general-`g` twins of `RMT/General/Trace.lean:325`'s `norm_stieltjesC_sub_le`
   and of the Cauchy estimate `RMT/General/Stability.lean:256`'s
   `norm_stieltjes2C_sub_mCDeriv_le`, proved directly from the Sherman-Morrison identity
   `R4C.resolvC_sub_vecMulVec` (`RMT/General/Companion.lean:176`) instead of through
   `gram`/`updateRow`.
3. Namespace `SpikedModel`: the five general theorems, each the twin of the matching
   Gaussian one of `RMT/R3minus.lean:278` to `:355` with `hG : m.GaussianNoise` replaced by
   `hν : NoiseLaw ν`, `hG : m.GeneralNoise ν`, and `hn2 : ∀ N, 2 ≤ n N` dropped.
-/

open Filter Topology MeasureTheory ProbabilityTheory
open scoped Matrix ENNReal

namespace StackedSVD
namespace GenRMT
namespace R3minus

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

/-! ### The rank-one downdate, at a general `W` and a general `g`

The Gaussian file routes R1b through R0's block law, which needs `2 ≤ n N`. Here we downdate
`GenRMT.gram (m.Z N ω)` directly by `Matrix.vecMulVec (m.gvec N ω) (m.gvec N ω)`
(`SpikedModel.W0_eq_gram_sub`, `RMT/General/Companion.lean:391`), so the two lemmas below are
stated at a general Hermitian `W` and a general vector `g`, not through `gram`/`updateRow`. -/

/-- **The rank-one downdate barely moves the trace, at every point of the upper half plane.**
General-`W`, general-`g` twin of `norm_stieltjesC_sub_le` (`RMT/General/Trace.lean:325`),
proved directly from the Sherman-Morrison identity `R4C.resolvC_sub_vecMulVec`
(`RMT/General/Companion.lean:176`) instead of through `gram`/`updateRow`. Stated at every `ζ`
in the upper half plane, not only at a fixed `z`, because the Cauchy estimate below applies it
on a whole circle. -/
theorem norm_stieltjesC_sub_vecMulVec_le {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ}
    (hW : W.IsHermitian) (hD : 0 < D) (g : Fin D → ℝ) {ζ : ℂ} (hζ : 0 < ζ.im) :
    ‖R4C.stieltjesC (W - Matrix.vecMulVec g g) ζ - R4C.stieltjesC W ζ‖
      ≤ 1 / ((D : ℝ) * ζ.im) := by
  have hζ' : ζ.im ≠ 0 := hζ.ne'
  have hsymm : R4C.cvec g ⬝ᵥ (R4C.resolvC W ζ *ᵥ (R4C.resolvC W ζ *ᵥ R4C.cvec g))
      = (R4C.resolvC W ζ *ᵥ R4C.cvec g) ⬝ᵥ (R4C.resolvC W ζ *ᵥ R4C.cvec g) := by
    rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose,
      ResolvDeriv.transpose_resolvC hW hζ']
  have hww : (R4C.resolvC W ζ *ᵥ R4C.cvec g) ⬝ᵥ (R4C.resolvC W ζ *ᵥ R4C.cvec g)
      = R4C.qform2C W ζ g := by
    rw [← hsymm, Matrix.mulVec_mulVec]
    rfl
  have htrvv : (Matrix.vecMulVec (R4C.resolvC W ζ *ᵥ R4C.cvec g)
      (R4C.resolvC W ζ *ᵥ R4C.cvec g)).trace = R4C.qform2C W ζ g := by
    rw [← hww]
    simp [Matrix.trace, Matrix.diag, Matrix.vecMulVec_apply, dotProduct]
  have htrace : (R4C.resolvC (W - Matrix.vecMulVec g g) ζ).trace - (R4C.resolvC W ζ).trace
      = (1 - R4C.qformC W ζ g)⁻¹ * R4C.qform2C W ζ g := by
    rw [R4C.resolvC_sub_vecMulVec hW hζ g, Matrix.trace_add, Matrix.trace_smul, smul_eq_mul,
      htrvv]
    ring
  have hne : (0 : ℝ) < ‖1 - R4C.qformC W ζ g‖ :=
    norm_pos_iff.mpr (sub_ne_zero.mpr (Ne.symm (R4C.qformC_ne_one hW hζ g)))
  have hbnd : ζ.im * ‖R4C.qform2C W ζ g‖ ≤ ‖1 - R4C.qformC W ζ g‖ :=
    im_mul_norm_qform2C_le hW hζ g
  have htrbound : ‖(R4C.resolvC (W - Matrix.vecMulVec g g) ζ).trace
      - (R4C.resolvC W ζ).trace‖ ≤ 1 / ζ.im := by
    rw [htrace, norm_mul, norm_inv, inv_mul_eq_div, div_le_div_iff₀ hne hζ, mul_comm]
    linarith
  have hDR : (0 : ℝ) < D := by exact_mod_cast hD
  have hstep2 : R4C.stieltjesC (W - Matrix.vecMulVec g g) ζ - R4C.stieltjesC W ζ
      = (D : ℂ)⁻¹ * ((R4C.resolvC (W - Matrix.vecMulVec g g) ζ).trace
          - (R4C.resolvC W ζ).trace) := by
    simp only [R4C.stieltjesC]
    ring
  rw [hstep2, norm_mul, norm_inv, Complex.norm_natCast]
  calc ((D : ℝ))⁻¹ * ‖(R4C.resolvC (W - Matrix.vecMulVec g g) ζ).trace
        - (R4C.resolvC W ζ).trace‖
      ≤ ((D : ℝ))⁻¹ * (1 / ζ.im) := mul_le_mul_of_nonneg_left htrbound (by positivity)
    _ = 1 / ((D : ℝ) * ζ.im) := by field_simp


/-- **The rank-one downdate barely moves `stieltjes2C`, at `z`.** General-`W`, general-`g`
twin of the Cauchy estimate `norm_stieltjes2C_sub_mCDeriv_le` (`RMT/General/Stability.lean:256`),
with `R4C.stieltjesC W` in place of `MP.mC c`: the difference function is holomorphic on the
whole upper half plane, `norm_stieltjesC_sub_vecMulVec_le` bounds it by `1/(D η)` at every
point, so on the circle of radius `z.im / 2` around `z` it is at most `2/(D z.im)`, and
Cauchy's estimate turns that uniform bound into `4/(D z.im²)` on the derivative gap at the
center. -/
theorem norm_stieltjes2C_sub_vecMulVec_le {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ}
    (hW : W.IsHermitian) {z : ℂ} (hz : 0 < z.im) (g : Fin D → ℝ) (hD : 0 < D) :
    ‖R4C.stieltjes2C (W - Matrix.vecMulVec g g) z - R4C.stieltjes2C W z‖
      ≤ 4 / ((D : ℝ) * z.im ^ 2) := by
  have hDR : (0 : ℝ) < D := by exact_mod_cast hD
  have hWvv : (W - Matrix.vecMulVec g g).IsHermitian := by
    have hvv : (Matrix.vecMulVec g g : Matrix (Fin D) (Fin D) ℝ).IsHermitian := by
      ext i j
      simp [Matrix.conjTranspose_apply, Matrix.vecMulVec_apply, mul_comm]
    exact hW.sub hvv
  have hrpos : (0 : ℝ) < z.im / 2 := by linarith
  have hball : ∀ ζ ∈ Metric.closedBall z (z.im / 2), 0 < ζ.im := by
    intro ζ hζ
    have := im_ge_of_mem_closedBall (z := z) hζ
    linarith
  have hderivfun : ∀ ζ : ℂ, 0 < ζ.im →
      HasDerivAt (fun w => R4C.stieltjesC (W - Matrix.vecMulVec g g) w - R4C.stieltjesC W w)
        (R4C.stieltjes2C (W - Matrix.vecMulVec g g) ζ - R4C.stieltjes2C W ζ) ζ :=
    fun ζ hζ => (hasDerivAt_stieltjesC hWvv hζ).sub (hasDerivAt_stieltjesC hW hζ)
  have hdiffOn : DifferentiableOn ℂ
      (fun w => R4C.stieltjesC (W - Matrix.vecMulVec g g) w - R4C.stieltjesC W w)
      (Metric.closedBall z (z.im / 2)) := fun ζ hζ =>
    ((hderivfun ζ (hball ζ hζ)).differentiableAt).differentiableWithinAt
  have hdc : DiffContOnCl ℂ
      (fun w => R4C.stieltjesC (W - Matrix.vecMulVec g g) w - R4C.stieltjesC W w)
      (Metric.ball z (z.im / 2)) := by
    constructor
    · exact hdiffOn.mono Metric.ball_subset_closedBall
    · rw [closure_ball z (ne_of_gt hrpos)]
      exact hdiffOn.continuousOn
  have hderiv : deriv
      (fun w => R4C.stieltjesC (W - Matrix.vecMulVec g g) w - R4C.stieltjesC W w) z
      = R4C.stieltjes2C (W - Matrix.vecMulVec g g) z - R4C.stieltjes2C W z :=
    (hderivfun z hz).deriv
  have hC : ∀ ζ ∈ Metric.sphere z (z.im / 2),
      ‖R4C.stieltjesC (W - Matrix.vecMulVec g g) ζ - R4C.stieltjesC W ζ‖
        ≤ 2 / ((D : ℝ) * z.im) := by
    intro ζ hζmem
    have hζball : ζ ∈ Metric.closedBall z (z.im / 2) := Metric.sphere_subset_closedBall hζmem
    have hζim : 0 < ζ.im := hball ζ hζball
    have hζlarge : z.im / 2 ≤ ζ.im := by
      have := im_ge_of_mem_closedBall (z := z) hζball
      linarith
    have h1 := norm_stieltjesC_sub_vecMulVec_le hW hD g hζim
    calc ‖R4C.stieltjesC (W - Matrix.vecMulVec g g) ζ - R4C.stieltjesC W ζ‖
        ≤ 1 / ((D : ℝ) * ζ.im) := h1
      _ ≤ 2 / ((D : ℝ) * z.im) := by
          rw [div_le_div_iff₀ (by positivity) (by positivity)]
          nlinarith [mul_le_mul_of_nonneg_left hζlarge hDR.le]
  rw [← hderiv]
  have hbound := Complex.norm_deriv_le_of_forall_mem_sphere_norm_le hrpos hdc hC
  have heq : (2 / ((D : ℝ) * z.im)) / (z.im / 2) = 4 / ((D : ℝ) * z.im ^ 2) := by
    have hzim : z.im ≠ 0 := hz.ne'
    have hDne : (D : ℝ) ≠ 0 := hDR.ne'
    field_simp
    ring
  rwa [heq] at hbound

end GenRMT

/-! ### The model: `W₀`, the Gram matrix, and the subcritical `lamMax` field

The three deterministic helper lemmas below carry no probability and no noise law; they are
the same content as `RMT/R3minus.lean:237` to `:269`, reproved here because `RMT/General/`
may not import that file (choice 8 of the plan note). -/

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- `ω ↦ X N ω` is measurable, entrywise, from `m.hZ`. Same route as `measurable_Eperp` of
`RMT/R5.lean`. -/
private theorem measurable_X (m : SpikedModel μ n d) (N : ℕ) : Measurable (m.X N) := by
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
private theorem measurableSet_gramLamMax_le (m : SpikedModel μ n d) (N : ℕ) (r : ℝ) :
    MeasurableSet {ω | gramLamMax (m.X N ω) ≤ r} := by
  have h : {ω | gramLamMax (m.X N ω) ≤ r}
      = (fun ω => gramLamMax (m.X N ω)) ⁻¹' Set.Iic r := rfl
  rw [h]
  exact (measurable_gramLamMax.comp (m.measurable_X N)) measurableSet_Iic

/-- The rank-one part only raises the top eigenvalue: `λ_max(W₀) ≤ λ_max(Xᵀ X)`. This is
`R4.lamMax_le_lamMax_vecMulVec` transported along R0's `gram_eq`. -/
private theorem lamMax_W0_le_gramLamMax (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) :
    lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ gramLamMax (m.X N ω) := by
  have hgram := m.gram_eq N ω
  have hH : (m.W0 N ω + Matrix.vecMulVec (m.qvec N ω) (m.qvec N ω)).IsHermitian := by
    rw [← hgram]; exact isHermitian_transpose_mul_self _
  calc lamMax (m.W0 N ω) (m.isHermitian_W0 N ω)
      ≤ lamMax (m.W0 N ω + Matrix.vecMulVec (m.qvec N ω) (m.qvec N ω)) hH :=
        R4.lamMax_le_lamMax_vecMulVec (m.isHermitian_W0 N ω) hH
    _ = gramLamMax (m.X N ω) :=
        (lamMax_congr hgram (isHermitian_transpose_mul_self (m.X N ω)) hH).symm

/-- A deterministic (`ω`-independent), nonnegative sequence tending to `0` classically tends
to `0` in probability. Used once below, to fold the `4/(d N z.im²)` error of the rank-one
downdate into the general-law limit of unit G3. Public since F37 (2026-09-09), the canonical
copy for `RMT/General/` (the private twin of `FormsLimits.lean` is dropped). -/
theorem tendstoInProb_of_tendsto_zero {a : ℕ → ℝ} (ha : ∀ N, 0 ≤ a N)
    (hatend : Tendsto a atTop (𝓝 0)) :
    TendstoInProb μ (fun N (_ : Ω N) => a N) 0 := by
  intro δ hδ
  have hev : ∀ᶠ N in atTop, a N < δ := hatend.eventually (gt_mem_nhds hδ)
  have heq : (fun _ : ℕ => (0 : ℝ≥0∞))
      =ᶠ[atTop] (fun N => μ N {ω : Ω N | δ ≤ |a N - 0|}) := by
    filter_upwards [hev] with N hN
    have hset : {ω : Ω N | δ ≤ |a N - 0|} = (∅ : Set (Ω N)) := by
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false, sub_zero,
        abs_of_nonneg (ha N)]
      exact not_le.mpr hN
    rw [hset, measure_empty]
  exact (tendsto_congr' heq).mp tendsto_const_nhds

variable [∀ N, IsProbabilityMeasure (μ N)] {m : SpikedModel μ n d} {c : ℝ}

omit [∀ N, IsProbabilityMeasure (μ N)] in
/-- **The general-law twin of R1b on the model's `W₀`.** Feeds unit G3's
`GenRMT.tendstoInProb_stieltjes2C_general` at `Y := m.Z`, `pN := n`, `dN := d`, then bridges
`GenRMT.gram (m.Z N ω)` to `m.W0 N ω` by the rank-one downdate of unit G2
(`SpikedModel.W0_eq_gram_sub`) through the Cauchy estimate
`GenRMT.norm_stieltjes2C_sub_vecMulVec_le`, whose error `4/(d N · z.im²)` is a deterministic
null sequence (`hreg.2.1` sends `d N → ∞`). No `2 ≤ n N` is needed: the block of R0 is never
formed, which is the real gain of the general route. -/
theorem tendstoInProb_stieltjes2C_W0_general (hc : 0 < c) (hreg : m.Regime c)
    {ν : Measure ℝ} (hν : NoiseLaw ν) (hG : m.GeneralNoise ν) {z : ℂ} (hz : 0 < z.im) :
    TendstoInProb μ (fun N ω => ‖R4C.stieltjes2C (m.W0 N ω) z - MP.mCDeriv c z‖) 0 := by
  have hgen := GenRMT.tendstoInProb_stieltjes2C_general (c := c) (z := z) (pN := n) (dN := d)
    hν hc hz hreg.2.1 m.hn hreg.2.2 m.Z hG
  have hdetTend : Tendsto (fun N => (4 : ℝ) / ((d N : ℝ) * z.im ^ 2)) atTop (𝓝 0) := by
    have hinv : Tendsto (fun N => ((d N : ℝ))⁻¹) atTop (𝓝 0) :=
      (tendsto_natCast_atTop_atTop.comp hreg.2.1).inv_tendsto_atTop
    have hmul := hinv.const_mul (4 / z.im ^ 2)
    rw [mul_zero] at hmul
    refine hmul.congr fun N => ?_
    rcases eq_or_ne ((d N : ℝ)) 0 with h0 | h0
    · rw [h0]; simp
    · field_simp
  have hdetnn : ∀ N, (0 : ℝ) ≤ (4 : ℝ) / ((d N : ℝ) * z.im ^ 2) := fun N => by positivity
  have hdetTIP : TendstoInProb μ (fun N (_ : Ω N) => (4 : ℝ) / ((d N : ℝ) * z.im ^ 2)) 0 :=
    tendstoInProb_of_tendsto_zero hdetnn hdetTend
  have hsum : TendstoInProb μ
      (fun N ω => ‖R4C.stieltjes2C (GenRMT.gram (m.Z N ω)) z - MP.mCDeriv c z‖
        + (4 : ℝ) / ((d N : ℝ) * z.im ^ 2)) 0 := by
    simpa using hgen.add hdetTIP
  refine TendstoInProb.of_le (fun N => Filter.Eventually.of_forall fun ω => ?_) hsum
  rw [sub_zero, abs_of_nonneg (norm_nonneg _), m.W0_eq_gram_sub N ω]
  have hbnd := GenRMT.norm_stieltjes2C_sub_vecMulVec_le (GenRMT.gram_isHermitian (m.Z N ω))
    hz (m.gvec N ω) (m.hd N)
  have e : R4C.stieltjes2C (GenRMT.gram (m.Z N ω)
        - Matrix.vecMulVec (m.gvec N ω) (m.gvec N ω)) z - MP.mCDeriv c z
      = (R4C.stieltjes2C (GenRMT.gram (m.Z N ω)
            - Matrix.vecMulVec (m.gvec N ω) (m.gvec N ω)) z
          - R4C.stieltjes2C (GenRMT.gram (m.Z N ω)) z)
        + (R4C.stieltjes2C (GenRMT.gram (m.Z N ω)) z - MP.mCDeriv c z) := by ring
  rw [e]
  calc ‖(R4C.stieltjes2C (GenRMT.gram (m.Z N ω)
          - Matrix.vecMulVec (m.gvec N ω) (m.gvec N ω)) z
        - R4C.stieltjes2C (GenRMT.gram (m.Z N ω)) z)
        + (R4C.stieltjes2C (GenRMT.gram (m.Z N ω)) z - MP.mCDeriv c z)‖
      ≤ ‖R4C.stieltjes2C (GenRMT.gram (m.Z N ω)
            - Matrix.vecMulVec (m.gvec N ω) (m.gvec N ω)) z
          - R4C.stieltjes2C (GenRMT.gram (m.Z N ω)) z‖
        + ‖R4C.stieltjes2C (GenRMT.gram (m.Z N ω)) z - MP.mCDeriv c z‖ := norm_add_le _ _
    _ ≤ (4 : ℝ) / ((d N : ℝ) * z.im ^ 2)
        + ‖R4C.stieltjes2C (GenRMT.gram (m.Z N ω)) z - MP.mCDeriv c z‖ := add_le_add hbnd le_rfl
    _ = ‖R4C.stieltjes2C (GenRMT.gram (m.Z N ω)) z - MP.mCDeriv c z‖
        + (4 : ℝ) / ((d N : ℝ) * z.im ^ 2) := by ring


/-- The trace form of (H2) at a real `x > bulkEdge c`, for the model's `W₀`: step 3 at the
edge bound of `ResolventLimits` and at the general-law twin of R1b,
`tendstoInProb_stieltjes2C_W0_general`. -/
theorem tendstoInProb_trace2_W0_general (H : m.ResolventLimits c) (hc : 0 < c)
    (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν) (hG : m.GeneralNoise ν)
    {x : ℝ} (hx : bulkEdge c < x) :
    TendstoInProb μ (fun N ω =>
      (d N : ℝ)⁻¹ * (R4.resolv (m.W0 N ω) x * R4.resolv (m.W0 N ω) x).trace)
      (MP.mDeriv c x) := by
  refine GenRMT.R3minus.tendstoInProb_trace2_of_complex hc (fun N ω => m.W0 N ω)
    (fun N ω => m.isHermitian_W0 N ω) (fun N => m.hd N) (fun ε hε => ?_)
    (fun z hz => tendstoInProb_stieltjes2C_W0_general hc hreg hν hG hz) hx
  exact tendsto_measure_compl_zero
    (fun N => (m.measurableSet_lamMax_W0_le N _).nullMeasurableSet) (H.edge ε hε)

/-- **R3⁻ for the model, general law.** The top eigenvalue of `W₀` does not fall below
`bulkEdge c`. -/
theorem tendsto_measure_lamMax_W0_le_sub_general (H : m.ResolventLimits c) (hc : 0 < c)
    (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν) (hG : m.GeneralNoise ν)
    {δ : ℝ} (hδ : 0 < δ) :
    Tendsto (fun N => μ N {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ bulkEdge c - δ})
      atTop (𝓝 0) :=
  GenRMT.R3minus.tendsto_measure_lamMax_le_sub hc (fun N ω => m.W0 N ω)
    (fun N ω => m.isHermitian_W0 N ω)
    (fun _ hx => tendstoInProb_trace2_W0_general H hc hreg hν hG hx) hδ

/-- **R3⁻ for the Gram matrix, general law.** `λ_max(Xᵀ X) ≥ λ_max(W₀)`, so the Gram top
eigenvalue does not fall below `bulkEdge c` either. -/
theorem tendsto_measure_gramLamMax_le_sub_general (H : m.ResolventLimits c) (hc : 0 < c)
    (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν) (hG : m.GeneralNoise ν)
    {δ : ℝ} (hδ : 0 < δ) :
    Tendsto (fun N => μ N {ω | gramLamMax (m.X N ω) ≤ bulkEdge c - δ}) atTop (𝓝 0) :=
  tendsto_measure_zero_of_subset
    (fun N _ hω => le_trans (m.lamMax_W0_le_gramLamMax N _) hω)
    (tendsto_measure_lamMax_W0_le_sub_general H hc hreg hν hG hδ)

/-- **The subcritical `lamMax` field of `SingleTableLaw`, general law.** For `θ⁴ ≤ c` (the
equality case included) `rhoSq θ c = bulkEdge c`, and the two edges pin the Gram top
eigenvalue there: R5's `tendsto_measure_lamMax_le_of_subcritical` gives the upper bound, this
file's general-law R3⁻ the lower one. No `2 ≤ n N` is needed anywhere in the chain: the real
gain of the general route over the Gaussian one. -/
theorem lamMax_tendstoInProb_of_subcritical_general (H : m.ResolventLimits c) (hc : 0 < c)
    (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν) (hG : m.GeneralNoise ν)
    (hθ : m.θ ^ 4 ≤ c) :
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
  have hlo := tendsto_measure_gramLamMax_le_sub_general H hc hreg hν hG hε
  refine tendsto_measure_zero_of_subset (fun N ω hω => ?_)
    (tendsto_measure_zero_union hup hlo)
  have hω' : ε ≤ |gramLamMax (m.X N ω) - bulkEdge c| := hω
  rcases le_abs.mp hω' with h | h
  · exact Set.mem_union_left _ (show bulkEdge c + ε ≤ gramLamMax (m.X N ω) by linarith)
  · exact Set.mem_union_right _ (show gramLamMax (m.X N ω) ≤ bulkEdge c - ε by linarith)

end SpikedModel

end StackedSVD
