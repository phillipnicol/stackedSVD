/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.R4C
import StackedSVD.RMT.MP7
import StackedSVD.Prob.TendstoInProb

/-!
# Item T: transfer from a complex `z` to the real axis

Review note: `notes/archive/rmt_T.md` (restructured 2026-08-30; decisions D6, D7, D9; choices 21,
24, 32). Items R1 and R2 give limits at `z` with `Im z = η > 0`, where `‖G₀(z)‖ ≤ 1/η` holds for
every realization. Item R5 needs the same limits at a real `z > bulkEdge c`. This file buys that
passage with an explicit error `O(η)`, so no Vitali theorem and no Stieltjes inversion enter.

## Content

* **T1, T2** (`norm_cform_sub_cformC_le`, `norm_cform2_sub_cform2C_le`): the deterministic
  bounds `4η √(v⬝v) √(w⬝w) / ε²` and `16η √(v⬝v) √(w⬝w) / ε³` between the real bilinear form
  at `x` and the complex one at `x + iη`, under `lamMax W ≤ bulkEdge c + ε/2` and
  `bulkEdge c + ε ≤ x`.
* **T1t, T2t** (`norm_trace_sub_stieltjesC_le`, `norm_trace2_sub_stieltjes2C_le`): the same
  for the two traces, with the constants `4η/ε²` and `16η/ε³` (the weights are `1/d` and sum
  to `1`).
* **T3 is deleted** (choice 21). `MP.mC_ofReal`, `MP.quad_mC`, `MP.im_mC_pos`,
  `MP.mC_eq_root`, `MP.tendsto_mC` and `MP.tendsto_mCDeriv` of `RMT/MP7.lean` are its six
  statements, each proved and each at least as strong.
* **T4** (`tendstoInProb_cform_of_complex`, `tendstoInProb_cform2_of_complex`, and the primed
  complement forms): convergence in probability of the real form at `x` from convergence of
  the complex form at every `x + iη`, for random vector sequences.

## The edge as a parameter (task H9, 2026-08-30)

No step of this file uses a property of `bulkEdge c`. The edge enters through the two real
inequalities `lamMax W ≤ b + ε/2` and `b + ε ≤ x` only. Namespace `T.Gen` therefore carries
every statement below with `bulkEdge c` replaced by a real parameter `b`, and with no
hypothesis on `b`. The eight names outside `T.Gen` are the MP case `b = bulkEdge c`; they
keep their 2026-08-30 signature and meaning, so no caller changes. The heteroscedastic
chain (`RMT/Het/`) calls `T.Gen.*` with `b = MPhet.bHet c w` or `b = MPhet.bSF c w`. The
two random vector sequences are `p` and `q` in `T.Gen` (`a` and `b` outside), because `b`
now names the edge.

## Statement change (recorded in `notes/archive/agent_reports/proof_t.md`)

The note states the T4 events as `μ N {good} → 1`. A measure is only an outer measure on a
set that is not measurable, so `μ N {good} → 1` does **not** give `μ N {good}ᶜ → 0` by
itself. The two T4 theorems therefore carry `[∀ N, IsProbabilityMeasure (μ N)]` and one
`NullMeasurableSet` hypothesis per event. The primed forms
`tendstoInProb_cform_of_complex'` and `tendstoInProb_cform2_of_complex'` state the same
conclusion with the complement events and need neither.
-/

open Filter Topology MeasureTheory
open scoped Matrix ENNReal

namespace StackedSVD
namespace T

variable {d : ℕ} {W : Matrix (Fin d) (Fin d) ℝ} {b c ε x η : ℝ} {v w : Fin d → ℝ}

/-! ### A sum estimate shared by the four deterministic bounds

Every one of T1, T2, T1t and T2t is a sum `∑ a, f a * (K a : ℂ)` against `∑ a, g a * (K a : ℂ)`
with a termwise bound `‖f a - g a‖ ≤ C` and a Cauchy-Schwarz bound `∑ a, |K a| ≤ S`. -/

private theorem norm_sum_sub_sum_le {f g : Fin d → ℂ} {K : Fin d → ℝ} {C S : ℝ}
    (hC : ∀ a, ‖f a - g a‖ ≤ C) (hCnn : 0 ≤ C) (hS : ∑ a, |K a| ≤ S) :
    ‖(∑ a, f a * (K a : ℂ)) - ∑ a, g a * (K a : ℂ)‖ ≤ C * S := by
  rw [← Finset.sum_sub_distrib]
  refine (norm_sum_le _ _).trans ?_
  have hterm : ∀ a ∈ (Finset.univ : Finset (Fin d)),
      ‖f a * (K a : ℂ) - g a * (K a : ℂ)‖ ≤ C * |K a| := by
    intro a _
    rw [← sub_mul, norm_mul, Complex.norm_real, Real.norm_eq_abs]
    exact mul_le_mul_of_nonneg_right (hC a) (abs_nonneg _)
  refine (Finset.sum_le_sum hterm).trans ?_
  rw [← Finset.mul_sum]
  exact mul_le_mul_of_nonneg_left hS hCnn

/-! ### The gap and the two termwise scalar bounds -/

/-- The edge hypothesis and `b + ε ≤ x` give the gap `ε/2` at every eigenvalue. -/
private theorem gap_ge (hW : W.IsHermitian) (hlam : lamMax W hW ≤ b + ε / 2)
    (hx : b + ε ≤ x) (a : Fin d) : ε / 2 ≤ x - hW.eigenvalues a := by
  have h := R4.eigenvalues_le_lamMax hW a
  linarith

private theorem ofReal_eigen_sub (hW : W.IsHermitian) (a : Fin d) :
    ((hW.eigenvalues a - x : ℝ) : ℂ) = ((-(x - hW.eigenvalues a) : ℝ) : ℂ) := by
  push_cast
  ring

private theorem eigen_sub_z (hW : W.IsHermitian) (a : Fin d) :
    (hW.eigenvalues a : ℂ) - ((x : ℂ) + (η : ℂ) * Complex.I)
      = ((-(x - hW.eigenvalues a) : ℝ) : ℂ) - (η : ℂ) * Complex.I := by
  push_cast
  ring

/-- **T1, termwise.** `|(λ - x)⁻¹ - (λ - x - iη)⁻¹| ≤ 4η/ε²`. -/
private theorem scalar_bound₁ (hW : W.IsHermitian) (hε : 0 < ε) (hη : 0 < η)
    (hlam : lamMax W hW ≤ b + ε / 2) (hx : b + ε ≤ x) (a : Fin d) :
    ‖((hW.eigenvalues a - x : ℝ) : ℂ)⁻¹
        - ((hW.eigenvalues a : ℂ) - ((x : ℂ) + (η : ℂ) * Complex.I))⁻¹‖ ≤ 4 * η / ε ^ 2 := by
  have hεne : ε ≠ 0 := hε.ne'
  rw [ofReal_eigen_sub hW a, eigen_sub_z hW a]
  refine (R4C.norm_inv_sub_inv_le (a := x - hW.eigenvalues a) (η := η) (g := ε / 2) hη
    (by linarith) (gap_ge hW hlam hx a)).trans (le_of_eq ?_)
  field_simp
  ring

/-- **T2, termwise.** `|(λ - x)⁻² - (λ - x - iη)⁻²| ≤ 16η/ε³`. -/
private theorem scalar_bound₂ (hW : W.IsHermitian) (hε : 0 < ε) (hη : 0 < η)
    (hlam : lamMax W hW ≤ b + ε / 2) (hx : b + ε ≤ x) (a : Fin d) :
    ‖(((hW.eigenvalues a - x : ℝ) : ℂ)⁻¹) ^ 2
        - (((hW.eigenvalues a : ℂ) - ((x : ℂ) + (η : ℂ) * Complex.I))⁻¹) ^ 2‖
      ≤ 16 * η / ε ^ 3 := by
  have hεne : ε ≠ 0 := hε.ne'
  rw [ofReal_eigen_sub hW a, eigen_sub_z hW a]
  refine (R4C.norm_inv_sq_sub_inv_sq_le (a := x - hW.eigenvalues a) (η := η) (g := ε / 2) hη
    (by linarith) (gap_ge hW hlam hx a)).trans (le_of_eq ?_)
  field_simp
  ring

namespace Gen

/-! ### T1 and T2: the two bilinear forms -/

set_option linter.unusedVariables false in
/-- **T1**, constant `4`, with the norm factor of choice 24. Cross forms are included;
`w = v` gives the quadratic form, since `R4.cform_self` and `R4C.qformC` are the diagonal
cases. -/
theorem norm_cform_sub_cformC_le (hW : W.IsHermitian) (hε : 0 < ε) (hη : 0 < η)
    (hlam : lamMax W hW ≤ b + ε / 2) (hx : b + ε ≤ x) :
    ‖((R4.cform W x v w : ℝ) : ℂ) - R4C.cformC W ((x : ℂ) + (η : ℂ) * Complex.I) v w‖
      ≤ 4 * η * Real.sqrt (v ⬝ᵥ v) * Real.sqrt (w ⬝ᵥ w) / ε ^ 2 := by
  have hlt : lamMax W hW < x := by linarith
  have hzne : ((x : ℂ) + (η : ℂ) * Complex.I).im ≠ 0 := by simpa using hη.ne'
  have hreal : ((R4.cform W x v w : ℝ) : ℂ)
      = ∑ a, ((hW.eigenvalues a - x : ℝ) : ℂ) ⁻¹ *
          ((((R4.eigU hW)ᵀ *ᵥ v) a * ((R4.eigU hW)ᵀ *ᵥ w) a : ℝ) : ℂ) := by
    have h : R4.cform W x v w
        = ∑ a, (hW.eigenvalues a - x)⁻¹ *
            (((R4.eigU hW)ᵀ *ᵥ v) a * ((R4.eigU hW)ᵀ *ᵥ w) a) := by
      rw [R4.cform, R4.resolv_eq_conj hW hlt, R4C.dotProduct_conj_gen]
    rw [h, Complex.ofReal_sum]
    exact Finset.sum_congr rfl fun a _ => by push_cast; ring
  rw [hreal, R4C.cformC_eq_sum hW hzne]
  have hS : ∑ a, |((R4.eigU hW)ᵀ *ᵥ v) a * ((R4.eigU hW)ᵀ *ᵥ w) a|
      ≤ Real.sqrt (v ⬝ᵥ v) * Real.sqrt (w ⬝ᵥ w) := by
    refine le_trans (le_of_eq (Finset.sum_congr rfl fun a _ => abs_mul _ _)) ?_
    exact R4C.sum_abs_coords_le hW v w
  refine (norm_sum_sub_sum_le (f := fun a => ((hW.eigenvalues a - x : ℝ) : ℂ)⁻¹)
    (g := fun a => ((hW.eigenvalues a : ℂ) - ((x : ℂ) + (η : ℂ) * Complex.I))⁻¹)
    (K := fun a => ((R4.eigU hW)ᵀ *ᵥ v) a * ((R4.eigU hW)ᵀ *ᵥ w) a)
    (fun a => scalar_bound₁ hW hε hη hlam hx a) (by positivity) hS).trans (le_of_eq ?_)
  ring

set_option linter.unusedVariables false in
/-- **T2**, constant `16`, same factor. -/
theorem norm_cform2_sub_cform2C_le (hW : W.IsHermitian) (hε : 0 < ε) (hη : 0 < η)
    (hlam : lamMax W hW ≤ b + ε / 2) (hx : b + ε ≤ x) :
    ‖((R4.cform2 W x v w : ℝ) : ℂ) - R4C.cform2C W ((x : ℂ) + (η : ℂ) * Complex.I) v w‖
      ≤ 16 * η * Real.sqrt (v ⬝ᵥ v) * Real.sqrt (w ⬝ᵥ w) / ε ^ 3 := by
  have hlt : lamMax W hW < x := by linarith
  have hzne : ((x : ℂ) + (η : ℂ) * Complex.I).im ≠ 0 := by simpa using hη.ne'
  have hreal : ((R4.cform2 W x v w : ℝ) : ℂ)
      = ∑ a, (((hW.eigenvalues a - x : ℝ) : ℂ)⁻¹) ^ 2 *
          ((((R4.eigU hW)ᵀ *ᵥ v) a * ((R4.eigU hW)ᵀ *ᵥ w) a : ℝ) : ℂ) := by
    have h : R4.cform2 W x v w
        = ∑ a, ((hW.eigenvalues a - x)⁻¹) ^ 2 *
            (((R4.eigU hW)ᵀ *ᵥ v) a * ((R4.eigU hW)ᵀ *ᵥ w) a) := by
      rw [R4.cform2, R4.resolv_mul_resolv_eq_conj hW hlt, R4C.dotProduct_conj_gen]
    rw [h, Complex.ofReal_sum]
    exact Finset.sum_congr rfl fun a _ => by push_cast; ring
  rw [hreal, R4C.cform2C_eq_sum hW hzne]
  have hS : ∑ a, |((R4.eigU hW)ᵀ *ᵥ v) a * ((R4.eigU hW)ᵀ *ᵥ w) a|
      ≤ Real.sqrt (v ⬝ᵥ v) * Real.sqrt (w ⬝ᵥ w) := by
    refine le_trans (le_of_eq (Finset.sum_congr rfl fun a _ => abs_mul _ _)) ?_
    exact R4C.sum_abs_coords_le hW v w
  refine (norm_sum_sub_sum_le (f := fun a => (((hW.eigenvalues a - x : ℝ) : ℂ)⁻¹) ^ 2)
    (g := fun a => (((hW.eigenvalues a : ℂ) - ((x : ℂ) + (η : ℂ) * Complex.I))⁻¹) ^ 2)
    (K := fun a => ((R4.eigU hW)ᵀ *ᵥ v) a * ((R4.eigU hW)ᵀ *ᵥ w) a)
    (fun a => scalar_bound₂ hW hε hη hlam hx a) (by positivity) hS).trans (le_of_eq ?_)
  ring

/-! ### T1t and T2t: the two trace forms -/

private theorem trace_resolv_eq_sum (hW : W.IsHermitian) (hlt : lamMax W hW < x) :
    (R4.resolv W x).trace = ∑ a, (hW.eigenvalues a - x)⁻¹ := by
  rw [R4.resolv_eq_conj hW hlt, Matrix.trace_mul_comm, ← Matrix.mul_assoc,
    R4.transpose_eigU_mul, Matrix.one_mul, Matrix.trace_diagonal]

private theorem trace_resolv_sq_eq_sum (hW : W.IsHermitian) (hlt : lamMax W hW < x) :
    (R4.resolv W x * R4.resolv W x).trace = ∑ a, ((hW.eigenvalues a - x)⁻¹) ^ 2 := by
  rw [R4.resolv_mul_resolv_eq_conj hW hlt, Matrix.trace_mul_comm, ← Matrix.mul_assoc,
    R4.transpose_eigU_mul, Matrix.one_mul, Matrix.trace_diagonal]

private theorem sum_abs_inv_card (hd : 0 < d) : ∑ _a : Fin d, |((d : ℝ))⁻¹| ≤ 1 := by
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd
  rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul,
    abs_of_nonneg (by positivity : (0 : ℝ) ≤ ((d : ℝ))⁻¹), mul_inv_cancel₀ hdR.ne']

set_option linter.unusedVariables false in
/-- **T1t.** Trace form, same constant. -/
theorem norm_trace_sub_stieltjesC_le (hW : W.IsHermitian) (hε : 0 < ε)
    (hη : 0 < η) (hd : 0 < d) (hlam : lamMax W hW ≤ b + ε / 2)
    (hx : b + ε ≤ x) :
    ‖(((d : ℝ)⁻¹ * (R4.resolv W x).trace : ℝ) : ℂ)
        - R4C.stieltjesC W ((x : ℂ) + (η : ℂ) * Complex.I)‖ ≤ 4 * η / ε ^ 2 := by
  have hlt : lamMax W hW < x := by linarith
  have hzne : ((x : ℂ) + (η : ℂ) * Complex.I).im ≠ 0 := by simpa using hη.ne'
  have hlhs : (((d : ℝ)⁻¹ * (R4.resolv W x).trace : ℝ) : ℂ)
      = ∑ a, ((hW.eigenvalues a - x : ℝ) : ℂ)⁻¹ * ((((d : ℝ))⁻¹ : ℝ) : ℂ) := by
    rw [trace_resolv_eq_sum hW hlt, Complex.ofReal_mul, Complex.ofReal_sum, Finset.mul_sum]
    exact Finset.sum_congr rfl fun a _ => by push_cast; ring
  have hrhs : R4C.stieltjesC W ((x : ℂ) + (η : ℂ) * Complex.I)
      = ∑ a, ((hW.eigenvalues a : ℂ) - ((x : ℂ) + (η : ℂ) * Complex.I))⁻¹ *
          ((((d : ℝ))⁻¹ : ℝ) : ℂ) := by
    rw [R4C.stieltjesC, R4C.trace_resolvC hW hzne, Finset.mul_sum]
    exact Finset.sum_congr rfl fun a _ => by push_cast; ring
  rw [hlhs, hrhs]
  refine (norm_sum_sub_sum_le (f := fun a => ((hW.eigenvalues a - x : ℝ) : ℂ)⁻¹)
    (g := fun a => ((hW.eigenvalues a : ℂ) - ((x : ℂ) + (η : ℂ) * Complex.I))⁻¹)
    (K := fun _ => ((d : ℝ))⁻¹)
    (fun a => scalar_bound₁ hW hε hη hlam hx a) (by positivity)
    (sum_abs_inv_card hd)).trans (le_of_eq ?_)
  ring

set_option linter.unusedVariables false in
/-- **T2t.** Trace form at second order; item R3⁻ uses this one. -/
theorem norm_trace2_sub_stieltjes2C_le (hW : W.IsHermitian) (hε : 0 < ε)
    (hη : 0 < η) (hd : 0 < d) (hlam : lamMax W hW ≤ b + ε / 2)
    (hx : b + ε ≤ x) :
    ‖(((d : ℝ)⁻¹ * (R4.resolv W x * R4.resolv W x).trace : ℝ) : ℂ)
        - R4C.stieltjes2C W ((x : ℂ) + (η : ℂ) * Complex.I)‖ ≤ 16 * η / ε ^ 3 := by
  have hlt : lamMax W hW < x := by linarith
  have hzne : ((x : ℂ) + (η : ℂ) * Complex.I).im ≠ 0 := by simpa using hη.ne'
  have hlhs : (((d : ℝ)⁻¹ * (R4.resolv W x * R4.resolv W x).trace : ℝ) : ℂ)
      = ∑ a, (((hW.eigenvalues a - x : ℝ) : ℂ)⁻¹) ^ 2 * ((((d : ℝ))⁻¹ : ℝ) : ℂ) := by
    rw [trace_resolv_sq_eq_sum hW hlt, Complex.ofReal_mul, Complex.ofReal_sum, Finset.mul_sum]
    exact Finset.sum_congr rfl fun a _ => by push_cast; ring
  have hrhs : R4C.stieltjes2C W ((x : ℂ) + (η : ℂ) * Complex.I)
      = ∑ a, (((hW.eigenvalues a : ℂ) - ((x : ℂ) + (η : ℂ) * Complex.I))⁻¹) ^ 2 *
          ((((d : ℝ))⁻¹ : ℝ) : ℂ) := by
    rw [R4C.stieltjes2C, R4C.trace_resolvC_sq hW hzne, Finset.mul_sum]
    exact Finset.sum_congr rfl fun a _ => by push_cast; ring
  rw [hlhs, hrhs]
  refine (norm_sum_sub_sum_le (f := fun a => (((hW.eigenvalues a - x : ℝ) : ℂ)⁻¹) ^ 2)
    (g := fun a => (((hW.eigenvalues a : ℂ) - ((x : ℂ) + (η : ℂ) * Complex.I))⁻¹) ^ 2)
    (K := fun _ => ((d : ℝ))⁻¹)
    (fun a => scalar_bound₂ hW hε hη hlam hx a) (by positivity)
    (sum_abs_inv_card hd)).trans (le_of_eq ?_)
  ring

end Gen

/-! ### The MP case of T1, T2, T1t and T2t: the old names at `b = bulkEdge c` -/

set_option linter.unusedVariables false in
/-- **T1** at the MP edge `b = bulkEdge c`. -/
theorem norm_cform_sub_cformC_le (hW : W.IsHermitian) (hc : 0 < c) (hε : 0 < ε) (hη : 0 < η)
    (hlam : lamMax W hW ≤ bulkEdge c + ε / 2) (hx : bulkEdge c + ε ≤ x) :
    ‖((R4.cform W x v w : ℝ) : ℂ) - R4C.cformC W ((x : ℂ) + (η : ℂ) * Complex.I) v w‖
      ≤ 4 * η * Real.sqrt (v ⬝ᵥ v) * Real.sqrt (w ⬝ᵥ w) / ε ^ 2 :=
  Gen.norm_cform_sub_cformC_le hW hε hη hlam hx

set_option linter.unusedVariables false in
/-- **T2** at the MP edge `b = bulkEdge c`. -/
theorem norm_cform2_sub_cform2C_le (hW : W.IsHermitian) (hc : 0 < c) (hε : 0 < ε) (hη : 0 < η)
    (hlam : lamMax W hW ≤ bulkEdge c + ε / 2) (hx : bulkEdge c + ε ≤ x) :
    ‖((R4.cform2 W x v w : ℝ) : ℂ) - R4C.cform2C W ((x : ℂ) + (η : ℂ) * Complex.I) v w‖
      ≤ 16 * η * Real.sqrt (v ⬝ᵥ v) * Real.sqrt (w ⬝ᵥ w) / ε ^ 3 :=
  Gen.norm_cform2_sub_cform2C_le hW hε hη hlam hx

set_option linter.unusedVariables false in
/-- **T1t** at the MP edge `b = bulkEdge c`. -/
theorem norm_trace_sub_stieltjesC_le (hW : W.IsHermitian) (hc : 0 < c) (hε : 0 < ε)
    (hη : 0 < η) (hd : 0 < d) (hlam : lamMax W hW ≤ bulkEdge c + ε / 2)
    (hx : bulkEdge c + ε ≤ x) :
    ‖(((d : ℝ)⁻¹ * (R4.resolv W x).trace : ℝ) : ℂ)
        - R4C.stieltjesC W ((x : ℂ) + (η : ℂ) * Complex.I)‖ ≤ 4 * η / ε ^ 2 :=
  Gen.norm_trace_sub_stieltjesC_le hW hε hη hd hlam hx

set_option linter.unusedVariables false in
/-- **T2t** at the MP edge `b = bulkEdge c`. -/
theorem norm_trace2_sub_stieltjes2C_le (hW : W.IsHermitian) (hc : 0 < c) (hε : 0 < ε)
    (hη : 0 < η) (hd : 0 < d) (hlam : lamMax W hW ≤ bulkEdge c + ε / 2)
    (hx : bulkEdge c + ε ≤ x) :
    ‖(((d : ℝ)⁻¹ * (R4.resolv W x * R4.resolv W x).trace : ℝ) : ℂ)
        - R4C.stieltjes2C W ((x : ℂ) + (η : ℂ) * Complex.I)‖ ≤ 16 * η / ε ^ 3 :=
  Gen.norm_trace2_sub_stieltjes2C_le hW hε hη hd hlam hx

/-! ### T3 is deleted (choice 21, 2026-08-30)

`MP.mC_ofReal`, `MP.quad_mC`, `MP.im_mC_pos`, `MP.mC_eq_root`, `MP.tendsto_mC` and
`MP.tendsto_mCDeriv` of `RMT/MP7.lean` are the six statements of the 2026-08-29 draft, each
proved and each at least as strong. Item T states nothing about `mC` and cites those names. -/

/-! ### T4: the transfer theorem, two forms, non-uniform -/

section Transfer

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {dN : ℕ → ℕ}

namespace Gen

/-- The choice of `η`: small enough for the T1 (or T2) error and for the scalar limit. -/
private theorem exists_eta {L : ℂ → ℂ} {ℓ x r : ℝ} (hr : 0 < r) {δ : ℝ} (hδ : 0 < δ)
    (hL : Tendsto (fun η : ℝ => L ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0) (𝓝 (ℓ : ℂ))) :
    ∃ η : ℝ, 0 < η ∧ η < r ∧ ‖L ((x : ℂ) + (η : ℂ) * Complex.I) - (ℓ : ℂ)‖ < δ := by
  have hball : Metric.ball (ℓ : ℂ) δ ∈ 𝓝 (ℓ : ℂ) := Metric.ball_mem_nhds _ hδ
  have hLnear : ∀ᶠ η : ℝ in 𝓝[>] (0 : ℝ),
      ‖L ((x : ℂ) + (η : ℂ) * Complex.I) - (ℓ : ℂ)‖ < δ := by
    filter_upwards [hL hball] with η hη
    simpa [Metric.mem_ball, dist_eq_norm] using hη
  have hsmall : ∀ᶠ η : ℝ in 𝓝[>] (0 : ℝ), η < r :=
    Filter.Eventually.filter_mono nhdsWithin_le_nhds
      (Filter.eventually_of_mem (Iio_mem_nhds hr) fun _ hy => hy)
  have hev : ∀ᶠ η : ℝ in 𝓝[>] (0 : ℝ),
      0 < η ∧ η < r ∧ ‖L ((x : ℂ) + (η : ℂ) * Complex.I) - (ℓ : ℂ)‖ < δ := by
    filter_upwards [self_mem_nhdsWithin, hsmall, hLnear] with η h1 h2 h3
    exact ⟨h1, h2, h3⟩
  obtain ⟨η, h1, h2, h3⟩ := hev.exists
  exact ⟨η, h1, h2, h3⟩

private theorem sqrt_le_two {y : ℝ} (hy : y ≤ 4) : Real.sqrt y ≤ 2 := by
  have h4 : Real.sqrt 4 = 2 := by
    rw [show (4 : ℝ) = 2 ^ 2 by norm_num]
    exact Real.sqrt_sq (by norm_num)
  calc Real.sqrt y ≤ Real.sqrt 4 := Real.sqrt_le_sqrt hy
    _ = 2 := h4

/-- **T4**, first order, complement form. `p` and `q` are random vector sequences, so the
`g`-forms of `ResolventLimits` go through the same theorem as the `v`-forms (mismatch M5).
`L` is the complex limit and `ℓ` its real axis value; instantiate `(L, ℓ)` by
`(MP.mC c, MP.m c x)` (`MP.tendsto_mC`) and by `(0, 0)`. -/
theorem tendstoInProb_cform_of_complex'
    (W₀ : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ)
    (hsymm : ∀ N ω, (W₀ N ω).IsHermitian)
    (p q : (N : ℕ) → Ω N → (Fin (dN N) → ℝ))
    (hedge : ∀ ε > 0, Tendsto
      (fun N => μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε}ᶜ) atTop (𝓝 0))
    (hnorm : Tendsto
      (fun N => μ N {ω | p N ω ⬝ᵥ p N ω ≤ 4 ∧ q N ω ⬝ᵥ q N ω ≤ 4}ᶜ) atTop (𝓝 0))
    (L : ℂ → ℂ) (ℓ : ℝ)
    (hcplx : ∀ z : ℂ, 0 < z.im →
      TendstoInProb μ (fun N ω => ‖R4C.cformC (W₀ N ω) z (p N ω) (q N ω) - L z‖) 0)
    {x : ℝ} (hx : b < x)
    (hL : Tendsto (fun η : ℝ => L ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0) (𝓝 (ℓ : ℂ))) :
    TendstoInProb μ (fun N ω => R4.cform (W₀ N ω) x (p N ω) (q N ω)) ℓ := by
  refine tendstoInProb_of_subset_union₃ fun δ hδ => ?_
  set ε : ℝ := x - b with hεdef
  have hε : 0 < ε := by rw [hεdef]; linarith
  have hxε : b + ε ≤ x := by rw [hεdef]; linarith
  obtain ⟨η, hη, hηsmall, hηL⟩ :=
    exists_eta (L := L) (ℓ := ℓ) (x := x) (r := δ * ε ^ 2 / 48) (by positivity)
      (by positivity : (0 : ℝ) < δ / 3) hL
  have hzim : (0 : ℝ) < ((x : ℂ) + (η : ℂ) * Complex.I).im := by simpa using hη
  refine ⟨fun N => {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε / 2}ᶜ,
    fun N => {ω | p N ω ⬝ᵥ p N ω ≤ 4 ∧ q N ω ⬝ᵥ q N ω ≤ 4}ᶜ,
    fun N => {ω | δ / 3 ≤ |‖R4C.cformC (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I)
      (p N ω) (q N ω) - L ((x : ℂ) + (η : ℂ) * Complex.I)‖ - 0|},
    ?_, hedge (ε / 2) (by positivity), hnorm,
    hcplx ((x : ℂ) + (η : ℂ) * Complex.I) hzim (δ / 3) (by positivity)⟩
  intro N ω hω
  by_contra hnot
  have hlam : lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε / 2 := by
    by_contra hcon
    exact hnot (Set.mem_union_left _ (Set.mem_union_left _ (Set.mem_compl hcon)))
  have hab : p N ω ⬝ᵥ p N ω ≤ 4 ∧ q N ω ⬝ᵥ q N ω ≤ 4 := by
    by_contra hcon
    exact hnot (Set.mem_union_left _ (Set.mem_union_right _ (Set.mem_compl hcon)))
  obtain ⟨haa, hbb⟩ := hab
  have hCf : ‖R4C.cformC (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I) (p N ω) (q N ω)
      - L ((x : ℂ) + (η : ℂ) * Complex.I)‖ < δ / 3 := by
    by_contra hcon
    have hcon' := not_lt.mp hcon
    refine hnot (Set.mem_union_right _ ?_)
    change δ / 3 ≤ |‖R4C.cformC (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I) (p N ω) (q N ω)
      - L ((x : ℂ) + (η : ℂ) * Complex.I)‖ - 0|
    rwa [sub_zero, abs_norm]
  have hT1 := norm_cform_sub_cformC_le (W := W₀ N ω) (v := p N ω) (w := q N ω)
    (hsymm N ω) hε hη hlam hxε
  have hna : Real.sqrt (p N ω ⬝ᵥ p N ω) ≤ 2 := sqrt_le_two haa
  have hnb : Real.sqrt (q N ω ⬝ᵥ q N ω) ≤ 2 := sqrt_le_two hbb
  have hstep : 4 * η * Real.sqrt (p N ω ⬝ᵥ p N ω) * Real.sqrt (q N ω ⬝ᵥ q N ω) / ε ^ 2
      < δ / 3 := by
    have hnum : 4 * η * Real.sqrt (p N ω ⬝ᵥ p N ω) * Real.sqrt (q N ω ⬝ᵥ q N ω)
        ≤ 16 * η := by
      have h2 : (0 : ℝ) ≤ Real.sqrt (q N ω ⬝ᵥ q N ω) := Real.sqrt_nonneg _
      have hprod : Real.sqrt (p N ω ⬝ᵥ p N ω) * Real.sqrt (q N ω ⬝ᵥ q N ω) ≤ 4 := by
        have := mul_le_mul hna hnb h2 (by norm_num : (0 : ℝ) ≤ 2)
        linarith
      have hmul := mul_le_mul_of_nonneg_left hprod (by positivity : (0 : ℝ) ≤ 4 * η)
      calc 4 * η * Real.sqrt (p N ω ⬝ᵥ p N ω) * Real.sqrt (q N ω ⬝ᵥ q N ω)
          = 4 * η * (Real.sqrt (p N ω ⬝ᵥ p N ω) * Real.sqrt (q N ω ⬝ᵥ q N ω)) := by ring
        _ ≤ 4 * η * 4 := hmul
        _ = 16 * η := by ring
    have hlt2 : 16 * η / ε ^ 2 < δ / 3 := by
      rw [div_lt_iff₀ (by positivity)]
      nlinarith [hηsmall, sq_nonneg ε]
    calc 4 * η * Real.sqrt (p N ω ⬝ᵥ p N ω) * Real.sqrt (q N ω ⬝ᵥ q N ω) / ε ^ 2
        ≤ 16 * η / ε ^ 2 := by
          exact div_le_div_of_nonneg_right hnum (by positivity)
      _ < δ / 3 := hlt2
  have habs : |R4.cform (W₀ N ω) x (p N ω) (q N ω) - ℓ|
      = ‖((R4.cform (W₀ N ω) x (p N ω) (q N ω) : ℝ) : ℂ) - ((ℓ : ℝ) : ℂ)‖ := by
    rw [← Complex.ofReal_sub, Complex.norm_real, Real.norm_eq_abs]
  have hδle : δ ≤ |R4.cform (W₀ N ω) x (p N ω) (q N ω) - ℓ| := hω
  rw [habs] at hδle
  have htri : ‖((R4.cform (W₀ N ω) x (p N ω) (q N ω) : ℝ) : ℂ) - ((ℓ : ℝ) : ℂ)‖
      ≤ ‖((R4.cform (W₀ N ω) x (p N ω) (q N ω) : ℝ) : ℂ)
          - R4C.cformC (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I) (p N ω) (q N ω)‖
        + (‖R4C.cformC (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I) (p N ω) (q N ω)
          - L ((x : ℂ) + (η : ℂ) * Complex.I)‖
        + ‖L ((x : ℂ) + (η : ℂ) * Complex.I) - ((ℓ : ℝ) : ℂ)‖) := by
    have e : ((R4.cform (W₀ N ω) x (p N ω) (q N ω) : ℝ) : ℂ) - ((ℓ : ℝ) : ℂ)
        = (((R4.cform (W₀ N ω) x (p N ω) (q N ω) : ℝ) : ℂ)
            - R4C.cformC (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I) (p N ω) (q N ω))
          + ((R4C.cformC (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I) (p N ω) (q N ω)
              - L ((x : ℂ) + (η : ℂ) * Complex.I))
            + (L ((x : ℂ) + (η : ℂ) * Complex.I) - ((ℓ : ℝ) : ℂ))) := by ring
    rw [e]
    exact (norm_add_le _ _).trans (add_le_add le_rfl (norm_add_le _ _))
  linarith [hT1, hstep, hCf, hηL, hδle, htri]

/-- **T4**, second order, complement form. Instantiate `(L, ℓ)` by
`(MP.mCDeriv c, MP.mDeriv c x)` (`MP.tendsto_mCDeriv`) and by `(0, 0)`. -/
theorem tendstoInProb_cform2_of_complex'
    (W₀ : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ)
    (hsymm : ∀ N ω, (W₀ N ω).IsHermitian)
    (p q : (N : ℕ) → Ω N → (Fin (dN N) → ℝ))
    (hedge : ∀ ε > 0, Tendsto
      (fun N => μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε}ᶜ) atTop (𝓝 0))
    (hnorm : Tendsto
      (fun N => μ N {ω | p N ω ⬝ᵥ p N ω ≤ 4 ∧ q N ω ⬝ᵥ q N ω ≤ 4}ᶜ) atTop (𝓝 0))
    (L : ℂ → ℂ) (ℓ : ℝ)
    (hcplx : ∀ z : ℂ, 0 < z.im →
      TendstoInProb μ (fun N ω => ‖R4C.cform2C (W₀ N ω) z (p N ω) (q N ω) - L z‖) 0)
    {x : ℝ} (hx : b < x)
    (hL : Tendsto (fun η : ℝ => L ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0) (𝓝 (ℓ : ℂ))) :
    TendstoInProb μ (fun N ω => R4.cform2 (W₀ N ω) x (p N ω) (q N ω)) ℓ := by
  refine tendstoInProb_of_subset_union₃ fun δ hδ => ?_
  set ε : ℝ := x - b with hεdef
  have hε : 0 < ε := by rw [hεdef]; linarith
  have hxε : b + ε ≤ x := by rw [hεdef]; linarith
  obtain ⟨η, hη, hηsmall, hηL⟩ :=
    exists_eta (L := L) (ℓ := ℓ) (x := x) (r := δ * ε ^ 3 / 192) (by positivity)
      (by positivity : (0 : ℝ) < δ / 3) hL
  have hzim : (0 : ℝ) < ((x : ℂ) + (η : ℂ) * Complex.I).im := by simpa using hη
  refine ⟨fun N => {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε / 2}ᶜ,
    fun N => {ω | p N ω ⬝ᵥ p N ω ≤ 4 ∧ q N ω ⬝ᵥ q N ω ≤ 4}ᶜ,
    fun N => {ω | δ / 3 ≤ |‖R4C.cform2C (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I)
      (p N ω) (q N ω) - L ((x : ℂ) + (η : ℂ) * Complex.I)‖ - 0|},
    ?_, hedge (ε / 2) (by positivity), hnorm,
    hcplx ((x : ℂ) + (η : ℂ) * Complex.I) hzim (δ / 3) (by positivity)⟩
  intro N ω hω
  by_contra hnot
  have hlam : lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε / 2 := by
    by_contra hcon
    exact hnot (Set.mem_union_left _ (Set.mem_union_left _ (Set.mem_compl hcon)))
  have hab : p N ω ⬝ᵥ p N ω ≤ 4 ∧ q N ω ⬝ᵥ q N ω ≤ 4 := by
    by_contra hcon
    exact hnot (Set.mem_union_left _ (Set.mem_union_right _ (Set.mem_compl hcon)))
  obtain ⟨haa, hbb⟩ := hab
  have hCf : ‖R4C.cform2C (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I) (p N ω) (q N ω)
      - L ((x : ℂ) + (η : ℂ) * Complex.I)‖ < δ / 3 := by
    by_contra hcon
    have hcon' := not_lt.mp hcon
    refine hnot (Set.mem_union_right _ ?_)
    change δ / 3 ≤ |‖R4C.cform2C (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I) (p N ω) (q N ω)
      - L ((x : ℂ) + (η : ℂ) * Complex.I)‖ - 0|
    rwa [sub_zero, abs_norm]
  have hT2 := norm_cform2_sub_cform2C_le (W := W₀ N ω) (v := p N ω) (w := q N ω)
    (hsymm N ω) hε hη hlam hxε
  have hna : Real.sqrt (p N ω ⬝ᵥ p N ω) ≤ 2 := sqrt_le_two haa
  have hnb : Real.sqrt (q N ω ⬝ᵥ q N ω) ≤ 2 := sqrt_le_two hbb
  have hstep : 16 * η * Real.sqrt (p N ω ⬝ᵥ p N ω) * Real.sqrt (q N ω ⬝ᵥ q N ω) / ε ^ 3
      < δ / 3 := by
    have hnum : 16 * η * Real.sqrt (p N ω ⬝ᵥ p N ω) * Real.sqrt (q N ω ⬝ᵥ q N ω)
        ≤ 64 * η := by
      have h2 : (0 : ℝ) ≤ Real.sqrt (q N ω ⬝ᵥ q N ω) := Real.sqrt_nonneg _
      have hprod : Real.sqrt (p N ω ⬝ᵥ p N ω) * Real.sqrt (q N ω ⬝ᵥ q N ω) ≤ 4 := by
        have := mul_le_mul hna hnb h2 (by norm_num : (0 : ℝ) ≤ 2)
        linarith
      have hmul := mul_le_mul_of_nonneg_left hprod (by positivity : (0 : ℝ) ≤ 16 * η)
      calc 16 * η * Real.sqrt (p N ω ⬝ᵥ p N ω) * Real.sqrt (q N ω ⬝ᵥ q N ω)
          = 16 * η * (Real.sqrt (p N ω ⬝ᵥ p N ω) * Real.sqrt (q N ω ⬝ᵥ q N ω)) := by ring
        _ ≤ 16 * η * 4 := hmul
        _ = 64 * η := by ring
    have hlt2 : 64 * η / ε ^ 3 < δ / 3 := by
      rw [div_lt_iff₀ (by positivity)]
      nlinarith [hηsmall, pow_pos hε 3]
    calc 16 * η * Real.sqrt (p N ω ⬝ᵥ p N ω) * Real.sqrt (q N ω ⬝ᵥ q N ω) / ε ^ 3
        ≤ 64 * η / ε ^ 3 := by
          exact div_le_div_of_nonneg_right hnum (by positivity)
      _ < δ / 3 := hlt2
  have habs : |R4.cform2 (W₀ N ω) x (p N ω) (q N ω) - ℓ|
      = ‖((R4.cform2 (W₀ N ω) x (p N ω) (q N ω) : ℝ) : ℂ) - ((ℓ : ℝ) : ℂ)‖ := by
    rw [← Complex.ofReal_sub, Complex.norm_real, Real.norm_eq_abs]
  have hδle : δ ≤ |R4.cform2 (W₀ N ω) x (p N ω) (q N ω) - ℓ| := hω
  rw [habs] at hδle
  have htri : ‖((R4.cform2 (W₀ N ω) x (p N ω) (q N ω) : ℝ) : ℂ) - ((ℓ : ℝ) : ℂ)‖
      ≤ ‖((R4.cform2 (W₀ N ω) x (p N ω) (q N ω) : ℝ) : ℂ)
          - R4C.cform2C (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I) (p N ω) (q N ω)‖
        + (‖R4C.cform2C (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I) (p N ω) (q N ω)
          - L ((x : ℂ) + (η : ℂ) * Complex.I)‖
        + ‖L ((x : ℂ) + (η : ℂ) * Complex.I) - ((ℓ : ℝ) : ℂ)‖) := by
    have e : ((R4.cform2 (W₀ N ω) x (p N ω) (q N ω) : ℝ) : ℂ) - ((ℓ : ℝ) : ℂ)
        = (((R4.cform2 (W₀ N ω) x (p N ω) (q N ω) : ℝ) : ℂ)
            - R4C.cform2C (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I) (p N ω) (q N ω))
          + ((R4C.cform2C (W₀ N ω) ((x : ℂ) + (η : ℂ) * Complex.I) (p N ω) (q N ω)
              - L ((x : ℂ) + (η : ℂ) * Complex.I))
            + (L ((x : ℂ) + (η : ℂ) * Complex.I) - ((ℓ : ℝ) : ℂ))) := by ring
    rw [e]
    exact (norm_add_le _ _).trans (add_le_add le_rfl (norm_add_le _ _))
  linarith [hT2, hstep, hCf, hηL, hδle, htri]

/-- **T4**, first order, in the shape of `notes/archive/rmt_T.md`: the two events carry probability
tending to `1`. The measure of a set that is not measurable is an outer measure, so the note's
`→ 1` needs a null measurable event and a probability measure to give `→ 0` on the
complement. -/
theorem tendstoInProb_cform_of_complex [∀ N, IsProbabilityMeasure (μ N)]
    (W₀ : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ)
    (hsymm : ∀ N ω, (W₀ N ω).IsHermitian)
    (p q : (N : ℕ) → Ω N → (Fin (dN N) → ℝ))
    (hedgeMeas : ∀ ε > 0, ∀ N,
      NullMeasurableSet {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε} (μ N))
    (hedge : ∀ ε > 0,
      Tendsto (fun N => μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε}) atTop (𝓝 1))
    (hnormMeas : ∀ N,
      NullMeasurableSet {ω | p N ω ⬝ᵥ p N ω ≤ 4 ∧ q N ω ⬝ᵥ q N ω ≤ 4} (μ N))
    (hnorm : Tendsto (fun N => μ N
      {ω | p N ω ⬝ᵥ p N ω ≤ 4 ∧ q N ω ⬝ᵥ q N ω ≤ 4}) atTop (𝓝 1))
    (L : ℂ → ℂ) (ℓ : ℝ)
    (hcplx : ∀ z : ℂ, 0 < z.im →
      TendstoInProb μ (fun N ω => ‖R4C.cformC (W₀ N ω) z (p N ω) (q N ω) - L z‖) 0)
    {x : ℝ} (hx : b < x)
    (hL : Tendsto (fun η : ℝ => L ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0) (𝓝 (ℓ : ℂ))) :
    TendstoInProb μ (fun N ω => R4.cform (W₀ N ω) x (p N ω) (q N ω)) ℓ :=
  tendstoInProb_cform_of_complex' W₀ hsymm p q
    (fun ε hε => tendsto_measure_compl_zero (hedgeMeas ε hε) (hedge ε hε))
    (tendsto_measure_compl_zero hnormMeas hnorm) L ℓ hcplx hx hL

/-- **T4**, second order, in the shape of `notes/archive/rmt_T.md`. -/
theorem tendstoInProb_cform2_of_complex [∀ N, IsProbabilityMeasure (μ N)]
    (W₀ : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ)
    (hsymm : ∀ N ω, (W₀ N ω).IsHermitian)
    (p q : (N : ℕ) → Ω N → (Fin (dN N) → ℝ))
    (hedgeMeas : ∀ ε > 0, ∀ N,
      NullMeasurableSet {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε} (μ N))
    (hedge : ∀ ε > 0,
      Tendsto (fun N => μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ b + ε}) atTop (𝓝 1))
    (hnormMeas : ∀ N,
      NullMeasurableSet {ω | p N ω ⬝ᵥ p N ω ≤ 4 ∧ q N ω ⬝ᵥ q N ω ≤ 4} (μ N))
    (hnorm : Tendsto (fun N => μ N
      {ω | p N ω ⬝ᵥ p N ω ≤ 4 ∧ q N ω ⬝ᵥ q N ω ≤ 4}) atTop (𝓝 1))
    (L : ℂ → ℂ) (ℓ : ℝ)
    (hcplx : ∀ z : ℂ, 0 < z.im →
      TendstoInProb μ (fun N ω => ‖R4C.cform2C (W₀ N ω) z (p N ω) (q N ω) - L z‖) 0)
    {x : ℝ} (hx : b < x)
    (hL : Tendsto (fun η : ℝ => L ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0) (𝓝 (ℓ : ℂ))) :
    TendstoInProb μ (fun N ω => R4.cform2 (W₀ N ω) x (p N ω) (q N ω)) ℓ :=
  tendstoInProb_cform2_of_complex' W₀ hsymm p q
    (fun ε hε => tendsto_measure_compl_zero (hedgeMeas ε hε) (hedge ε hε))
    (tendsto_measure_compl_zero hnormMeas hnorm) L ℓ hcplx hx hL

end Gen

/-! ### The MP case of T4: the old names at `b = bulkEdge c` -/

set_option linter.unusedVariables false in
/-- **T4, first order, complement form** at the MP edge `b = bulkEdge c`. -/
theorem tendstoInProb_cform_of_complex' (hc : 0 < c)
    (W₀ : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ)
    (hsymm : ∀ N ω, (W₀ N ω).IsHermitian)
    (a b : (N : ℕ) → Ω N → (Fin (dN N) → ℝ))
    (hedge : ∀ ε > 0, Tendsto
      (fun N => μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ bulkEdge c + ε}ᶜ) atTop (𝓝 0))
    (hnorm : Tendsto
      (fun N => μ N {ω | a N ω ⬝ᵥ a N ω ≤ 4 ∧ b N ω ⬝ᵥ b N ω ≤ 4}ᶜ) atTop (𝓝 0))
    (L : ℂ → ℂ) (ℓ : ℝ)
    (hcplx : ∀ z : ℂ, 0 < z.im →
      TendstoInProb μ (fun N ω => ‖R4C.cformC (W₀ N ω) z (a N ω) (b N ω) - L z‖) 0)
    {x : ℝ} (hx : bulkEdge c < x)
    (hL : Tendsto (fun η : ℝ => L ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0) (𝓝 (ℓ : ℂ))) :
    TendstoInProb μ (fun N ω => R4.cform (W₀ N ω) x (a N ω) (b N ω)) ℓ :=
  Gen.tendstoInProb_cform_of_complex' W₀ hsymm a b hedge hnorm L ℓ hcplx hx hL

set_option linter.unusedVariables false in
/-- **T4, second order, complement form** at the MP edge `b = bulkEdge c`. -/
theorem tendstoInProb_cform2_of_complex' (hc : 0 < c)
    (W₀ : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ)
    (hsymm : ∀ N ω, (W₀ N ω).IsHermitian)
    (a b : (N : ℕ) → Ω N → (Fin (dN N) → ℝ))
    (hedge : ∀ ε > 0, Tendsto
      (fun N => μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ bulkEdge c + ε}ᶜ) atTop (𝓝 0))
    (hnorm : Tendsto
      (fun N => μ N {ω | a N ω ⬝ᵥ a N ω ≤ 4 ∧ b N ω ⬝ᵥ b N ω ≤ 4}ᶜ) atTop (𝓝 0))
    (L : ℂ → ℂ) (ℓ : ℝ)
    (hcplx : ∀ z : ℂ, 0 < z.im →
      TendstoInProb μ (fun N ω => ‖R4C.cform2C (W₀ N ω) z (a N ω) (b N ω) - L z‖) 0)
    {x : ℝ} (hx : bulkEdge c < x)
    (hL : Tendsto (fun η : ℝ => L ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0) (𝓝 (ℓ : ℂ))) :
    TendstoInProb μ (fun N ω => R4.cform2 (W₀ N ω) x (a N ω) (b N ω)) ℓ :=
  Gen.tendstoInProb_cform2_of_complex' W₀ hsymm a b hedge hnorm L ℓ hcplx hx hL

set_option linter.unusedVariables false in
/-- **T4, first order, `→ 1` form** at the MP edge `b = bulkEdge c`. -/
theorem tendstoInProb_cform_of_complex [∀ N, IsProbabilityMeasure (μ N)] (hc : 0 < c)
    (W₀ : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ)
    (hsymm : ∀ N ω, (W₀ N ω).IsHermitian)
    (a b : (N : ℕ) → Ω N → (Fin (dN N) → ℝ))
    (hedgeMeas : ∀ ε > 0, ∀ N,
      NullMeasurableSet {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ bulkEdge c + ε} (μ N))
    (hedge : ∀ ε > 0,
      Tendsto (fun N => μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ bulkEdge c + ε}) atTop (𝓝 1))
    (hnormMeas : ∀ N,
      NullMeasurableSet {ω | a N ω ⬝ᵥ a N ω ≤ 4 ∧ b N ω ⬝ᵥ b N ω ≤ 4} (μ N))
    (hnorm : Tendsto (fun N => μ N
      {ω | a N ω ⬝ᵥ a N ω ≤ 4 ∧ b N ω ⬝ᵥ b N ω ≤ 4}) atTop (𝓝 1))
    (L : ℂ → ℂ) (ℓ : ℝ)
    (hcplx : ∀ z : ℂ, 0 < z.im →
      TendstoInProb μ (fun N ω => ‖R4C.cformC (W₀ N ω) z (a N ω) (b N ω) - L z‖) 0)
    {x : ℝ} (hx : bulkEdge c < x)
    (hL : Tendsto (fun η : ℝ => L ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0) (𝓝 (ℓ : ℂ))) :
    TendstoInProb μ (fun N ω => R4.cform (W₀ N ω) x (a N ω) (b N ω)) ℓ :=
  Gen.tendstoInProb_cform_of_complex W₀ hsymm a b hedgeMeas hedge hnormMeas hnorm L ℓ
    hcplx hx hL

set_option linter.unusedVariables false in
/-- **T4, second order, `→ 1` form** at the MP edge `b = bulkEdge c`. -/
theorem tendstoInProb_cform2_of_complex [∀ N, IsProbabilityMeasure (μ N)] (hc : 0 < c)
    (W₀ : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ)
    (hsymm : ∀ N ω, (W₀ N ω).IsHermitian)
    (a b : (N : ℕ) → Ω N → (Fin (dN N) → ℝ))
    (hedgeMeas : ∀ ε > 0, ∀ N,
      NullMeasurableSet {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ bulkEdge c + ε} (μ N))
    (hedge : ∀ ε > 0,
      Tendsto (fun N => μ N {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ bulkEdge c + ε}) atTop (𝓝 1))
    (hnormMeas : ∀ N,
      NullMeasurableSet {ω | a N ω ⬝ᵥ a N ω ≤ 4 ∧ b N ω ⬝ᵥ b N ω ≤ 4} (μ N))
    (hnorm : Tendsto (fun N => μ N
      {ω | a N ω ⬝ᵥ a N ω ≤ 4 ∧ b N ω ⬝ᵥ b N ω ≤ 4}) atTop (𝓝 1))
    (L : ℂ → ℂ) (ℓ : ℝ)
    (hcplx : ∀ z : ℂ, 0 < z.im →
      TendstoInProb μ (fun N ω => ‖R4C.cform2C (W₀ N ω) z (a N ω) (b N ω) - L z‖) 0)
    {x : ℝ} (hx : bulkEdge c < x)
    (hL : Tendsto (fun η : ℝ => L ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] 0) (𝓝 (ℓ : ℂ))) :
    TendstoInProb μ (fun N ω => R4.cform2 (W₀ N ω) x (a N ω) (b N ω)) ℓ :=
  Gen.tendstoInProb_cform2_of_complex W₀ hsymm a b hedgeMeas hedge hnormMeas hnorm L ℓ
    hcplx hx hL

end Transfer

end T
end StackedSVD
