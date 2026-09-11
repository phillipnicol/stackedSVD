/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Prob.TendstoInProb

/-!
# Convergence in probability for complex sequences

`StackedSVD.TendstoInProb` of `Defs.lean` is real valued. A complex limit is therefore
stated in norm form: `TendstoInProb μ (fun N ω => ‖f N ω - L‖) 0`. This file gives that
form the closure lemmas that `Prob/TendstoInProb.lean` gives the real form: sums,
differences, products, inverses, quotients, a constant, a congruence, a domination bound,
a triangle step through a random center, and the three maps back to a real limit (`norm`,
`re`, `im`).

Every statement here is a plain implication between limits in probability. Nothing needs a
probability measure, a measurable set, or a measurable function: each proof uses only
monotonicity and subadditivity of a measure, through the real API of
`Prob/TendstoInProb.lean`. The brief of unit F1a listed
`[∀ N, IsProbabilityMeasure (μ N)]` in the variable block; no proof uses it, so it is not
a hypothesis of any lemma here. A consumer that has the instance is unaffected.

`tendstoInProbC_iff` is the bridge that removes the `|· - 0|` wrapper. Use it whenever a
proof needs the event `{ω | ε ≤ ‖f N ω - L‖}` itself.

The one copy is `tendstoInProbC_trans`, from `RMT/R2.lean:1052`
(`tendstoInProb_norm_sub_trans`); `RMT.R2` is Gaussian only and this module must not import
it. The copy is renamed, not `private`, because the general Layer 2 route needs it.
-/

open MeasureTheory Filter Topology Set

namespace StackedSVD

section ProbC

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
variable {f g : ∀ N, Ω N → ℂ} {L M : ℂ}

/-! ### The bridge to the bare event -/

/-- A complex limit in probability, with the `|· - 0|` of the definition removed. Both
directions are one `simp` call: `|‖x‖ - 0| = ‖x‖`. -/
theorem tendstoInProbC_iff :
    TendstoInProb μ (fun N ω => ‖f N ω - L‖) 0 ↔
      ∀ ε > 0, Tendsto (fun N => μ N {ω | ε ≤ ‖f N ω - L‖}) atTop (𝓝 0) := by
  constructor <;> intro h ε hε <;> simpa using h ε hε

/-! ### Domination, congruence, constants -/

/-- Domination by a real sequence that tends to `0` in probability. -/
theorem tendstoInProbC_of_le {h : ∀ N, Ω N → ℝ} (hle : ∀ N ω, ‖f N ω - L‖ ≤ h N ω)
    (hh : TendstoInProb μ h 0) : TendstoInProb μ (fun N ω => ‖f N ω - L‖) 0 := by
  refine TendstoInProb.of_le (g := h) (fun N => Eventually.of_forall fun ω => ?_) hh
  rw [sub_zero, abs_of_nonneg (norm_nonneg _)]
  exact hle N ω

/-- A constant sequence. -/
theorem tendstoInProbC_const (μ : ∀ N, Measure (Ω N)) (L : ℂ) :
    TendstoInProb μ (fun _ (_ : Ω _) => ‖L - L‖) 0 := by
  simpa using TendstoInProb.const μ (0 : ℝ)

/-- Transfer along a pointwise identity. -/
theorem tendstoInProbC_congr (hfg : ∀ N ω, f N ω = g N ω)
    (hf : TendstoInProb μ (fun N ω => ‖f N ω - L‖) 0) :
    TendstoInProb μ (fun N ω => ‖g N ω - L‖) 0 := by
  refine TendstoInProb.congr hf fun N => Eventually.of_forall fun ω => ?_
  simp only [hfg]

/-- **Centering.** The gap to a random center plus the gap of the center. Copy of
`tendstoInProb_norm_sub_trans` of `RMT/R2.lean:1052`, renamed: `RMT.R2` is a Gaussian only
module and the general route must not import it. -/
theorem tendstoInProbC_trans {A B : ∀ N, Ω N → ℂ}
    (h1 : TendstoInProb μ (fun N ω => ‖A N ω - B N ω‖) 0)
    (h2 : TendstoInProb μ (fun N ω => ‖B N ω - L‖) 0) :
    TendstoInProb μ (fun N ω => ‖A N ω - L‖) 0 := by
  refine TendstoInProb.of_le (g := fun N ω => ‖A N ω - B N ω‖ + ‖B N ω - L‖) ?_ ?_
  · intro N
    filter_upwards with ω
    rw [sub_zero, abs_of_nonneg (norm_nonneg _)]
    calc ‖A N ω - L‖ = ‖A N ω - B N ω + (B N ω - L)‖ := by ring_nf
      _ ≤ ‖A N ω - B N ω‖ + ‖B N ω - L‖ := norm_add_le _ _
  · simpa using h1.add h2

/-! ### Algebra -/

/-- Sums. -/
theorem tendstoInProbC_add (hf : TendstoInProb μ (fun N ω => ‖f N ω - L‖) 0)
    (hg : TendstoInProb μ (fun N ω => ‖g N ω - M‖) 0) :
    TendstoInProb μ (fun N ω => ‖f N ω + g N ω - (L + M)‖) 0 := by
  refine tendstoInProbC_of_le (h := fun N ω => ‖f N ω - L‖ + ‖g N ω - M‖) (fun N ω => ?_) ?_
  · calc ‖f N ω + g N ω - (L + M)‖ = ‖f N ω - L + (g N ω - M)‖ := by ring_nf
      _ ≤ ‖f N ω - L‖ + ‖g N ω - M‖ := norm_add_le _ _
  · simpa using hf.add hg

/-- Negation. -/
theorem tendstoInProbC_neg (hf : TendstoInProb μ (fun N ω => ‖f N ω - L‖) 0) :
    TendstoInProb μ (fun N ω => ‖-f N ω - -L‖) 0 := by
  refine tendstoInProbC_of_le (h := fun N ω => ‖f N ω - L‖) (fun N ω => ?_) hf
  rw [neg_sub_neg, norm_sub_rev]

/-- Differences. -/
theorem tendstoInProbC_sub (hf : TendstoInProb μ (fun N ω => ‖f N ω - L‖) 0)
    (hg : TendstoInProb μ (fun N ω => ‖g N ω - M‖) 0) :
    TendstoInProb μ (fun N ω => ‖f N ω - g N ω - (L - M)‖) 0 := by
  refine tendstoInProbC_of_le (h := fun N ω => ‖f N ω - L‖ + ‖g N ω - M‖) (fun N ω => ?_) ?_
  · calc ‖f N ω - g N ω - (L - M)‖ = ‖f N ω - L - (g N ω - M)‖ := by ring_nf
      _ ≤ ‖f N ω - L‖ + ‖g N ω - M‖ := norm_sub_le _ _
  · simpa using hf.add hg

/-- A constant factor on the left. -/
theorem tendstoInProbC_const_mul (a : ℂ) (hf : TendstoInProb μ (fun N ω => ‖f N ω - L‖) 0) :
    TendstoInProb μ (fun N ω => ‖a * f N ω - a * L‖) 0 := by
  refine tendstoInProbC_of_le (h := fun N ω => ‖a‖ * ‖f N ω - L‖) (fun N ω => ?_) ?_
  · rw [← mul_sub, norm_mul]
  · simpa using TendstoInProb.const_mul ‖a‖ hf

/-- A constant factor on the right. -/
theorem tendstoInProbC_mul_const (a : ℂ) (hf : TendstoInProb μ (fun N ω => ‖f N ω - L‖) 0) :
    TendstoInProb μ (fun N ω => ‖f N ω * a - L * a‖) 0 := by
  refine tendstoInProbC_of_le (h := fun N ω => ‖f N ω - L‖ * ‖a‖) (fun N ω => ?_) ?_
  · rw [← sub_mul, norm_mul]
  · simpa using TendstoInProb.mul_const ‖a‖ hf

/-- Products. The bound is `‖fg - LM‖ ≤ ‖f - L‖ * (‖g - M‖ + ‖M‖) + ‖L‖ * ‖g - M‖`. -/
theorem tendstoInProbC_mul (hf : TendstoInProb μ (fun N ω => ‖f N ω - L‖) 0)
    (hg : TendstoInProb μ (fun N ω => ‖g N ω - M‖) 0) :
    TendstoInProb μ (fun N ω => ‖f N ω * g N ω - L * M‖) 0 := by
  refine tendstoInProbC_of_le
    (h := fun N ω => ‖f N ω - L‖ * (‖g N ω - M‖ + ‖M‖) + ‖L‖ * ‖g N ω - M‖)
    (fun N ω => ?_) ?_
  · have hsplit : f N ω * g N ω - L * M
        = (f N ω - L) * (g N ω - M) + (f N ω - L) * M + L * (g N ω - M) := by ring
    calc ‖f N ω * g N ω - L * M‖
        ≤ ‖(f N ω - L) * (g N ω - M) + (f N ω - L) * M‖ + ‖L * (g N ω - M)‖ := by
          rw [hsplit]; exact norm_add_le _ _
      _ ≤ ‖(f N ω - L) * (g N ω - M)‖ + ‖(f N ω - L) * M‖ + ‖L * (g N ω - M)‖ :=
          add_le_add (norm_add_le _ _) le_rfl
      _ = ‖f N ω - L‖ * (‖g N ω - M‖ + ‖M‖) + ‖L‖ * ‖g N ω - M‖ := by
          rw [norm_mul, norm_mul, norm_mul]; ring
  · simpa using (hf.mul (hg.add (TendstoInProb.const μ ‖M‖))).add
      (TendstoInProb.const_mul ‖L‖ hg)

/-- Inverses, for a nonzero limit. On the event `‖f - L‖ < ‖L‖ / 2` the value `f` is
nonzero and `‖f⁻¹ - L⁻¹‖ ≤ 2 * ‖f - L‖ / ‖L‖ ^ 2`; the complement has vanishing measure. -/
theorem tendstoInProbC_inv (hL : L ≠ 0) (hf : TendstoInProb μ (fun N ω => ‖f N ω - L‖) 0) :
    TendstoInProb μ (fun N ω => ‖(f N ω)⁻¹ - L⁻¹‖) 0 := by
  have hLpos : 0 < ‖L‖ := norm_pos_iff.mpr hL
  refine tendstoInProb_of_subset_union₂ fun δ hδ => ?_
  have hhalf : (0 : ℝ) < ‖L‖ / 2 := by linarith
  have hδL : (0 : ℝ) < δ * ‖L‖ ^ 2 / 2 :=
    div_pos (mul_pos hδ (pow_pos hLpos 2)) two_pos
  refine ⟨fun N => {ω | ‖L‖ / 2 ≤ ‖f N ω - L‖},
    fun N => {ω | δ * ‖L‖ ^ 2 / 2 ≤ ‖f N ω - L‖}, fun N ω hω => ?_,
    tendstoInProbC_iff.mp hf _ hhalf, tendstoInProbC_iff.mp hf _ hδL⟩
  have hω' : δ ≤ |‖(f N ω)⁻¹ - L⁻¹‖ - 0| := hω
  rw [sub_zero, abs_of_nonneg (norm_nonneg _)] at hω'
  by_cases hbig : ‖L‖ / 2 ≤ ‖f N ω - L‖
  · exact Set.mem_union_left _ hbig
  refine Set.mem_union_right _ ?_
  change δ * ‖L‖ ^ 2 / 2 ≤ ‖f N ω - L‖
  rw [not_le] at hbig
  -- `f N ω` is nonzero, with `‖L‖ / 2 < ‖f N ω‖`
  have hfnorm : ‖L‖ / 2 < ‖f N ω‖ := by
    have hrev := norm_sub_norm_le L (f N ω)
    rw [norm_sub_rev] at hrev
    linarith
  have hfne : f N ω ≠ 0 := norm_pos_iff.mp (by linarith)
  have hbound : ‖(f N ω)⁻¹ - L⁻¹‖ = ‖f N ω - L‖ / (‖f N ω‖ * ‖L‖) := by
    rw [inv_sub_inv hfne hL, norm_div, norm_mul, norm_sub_rev]
  have hden : (0 : ℝ) < ‖f N ω‖ * ‖L‖ := by
    exact mul_pos (by linarith) hLpos
  rw [hbound, le_div_iff₀ hden] at hω'
  have hprod : ‖L‖ / 2 * ‖L‖ ≤ ‖f N ω‖ * ‖L‖ :=
    mul_le_mul_of_nonneg_right hfnorm.le (norm_nonneg _)
  have hscale := mul_le_mul_of_nonneg_left hprod hδ.le
  have hring : δ * ‖L‖ ^ 2 / 2 = δ * (‖L‖ / 2 * ‖L‖) := by ring
  rw [hring]
  linarith [hscale, hω']

/-- Quotients, for a nonzero limit in the denominator. -/
theorem tendstoInProbC_div (hM : M ≠ 0) (hf : TendstoInProb μ (fun N ω => ‖f N ω - L‖) 0)
    (hg : TendstoInProb μ (fun N ω => ‖g N ω - M‖) 0) :
    TendstoInProb μ (fun N ω => ‖f N ω / g N ω - L / M‖) 0 := by
  simp only [div_eq_mul_inv]
  exact tendstoInProbC_mul hf (tendstoInProbC_inv hM hg)

/-! ### Between the real and the complex form -/

/-- A real limit in probability, read as a complex one. -/
theorem tendstoInProbC_ofReal {fr : ∀ N, Ω N → ℝ} {a : ℝ} (hf : TendstoInProb μ fr a) :
    TendstoInProb μ (fun N ω => ‖((fr N ω : ℂ)) - (a : ℂ)‖) 0 := by
  have habs : TendstoInProb μ (fun N ω => |fr N ω - a|) 0 := by
    intro ε hε
    simpa using hf ε hε
  refine tendstoInProbC_of_le (f := fun N ω => ((fr N ω : ℝ) : ℂ)) (L := (a : ℂ))
    (h := fun N ω => |fr N ω - a|) (fun N ω => le_of_eq ?_) habs
  rw [← Complex.ofReal_sub, Complex.norm_real, Real.norm_eq_abs]

/-- The norm of a complex sequence, as a real limit in probability. -/
theorem tendstoInProbC_norm (hf : TendstoInProb μ (fun N ω => ‖f N ω - L‖) 0) :
    TendstoInProb μ (fun N ω => ‖f N ω‖) ‖L‖ := by
  refine TendstoInProb.of_le (g := fun N ω => ‖f N ω - L‖)
    (fun N => Eventually.of_forall fun ω => ?_) hf
  exact abs_norm_sub_norm_le _ _

/-- The real part of a complex sequence, as a real limit in probability. -/
theorem tendstoInProbC_re (hf : TendstoInProb μ (fun N ω => ‖f N ω - L‖) 0) :
    TendstoInProb μ (fun N ω => (f N ω).re) L.re := by
  refine TendstoInProb.of_le (g := fun N ω => ‖f N ω - L‖)
    (fun N => Eventually.of_forall fun ω => ?_) hf
  rw [← Complex.sub_re]
  exact Complex.abs_re_le_norm _

/-- The imaginary part of a complex sequence, as a real limit in probability. -/
theorem tendstoInProbC_im (hf : TendstoInProb μ (fun N ω => ‖f N ω - L‖) 0) :
    TendstoInProb μ (fun N ω => (f N ω).im) L.im := by
  refine TendstoInProb.of_le (g := fun N ω => ‖f N ω - L‖)
    (fun N => Eventually.of_forall fun ω => ?_) hf
  rw [← Complex.sub_im]
  exact Complex.abs_im_le_norm _

end ProbC

end StackedSVD

