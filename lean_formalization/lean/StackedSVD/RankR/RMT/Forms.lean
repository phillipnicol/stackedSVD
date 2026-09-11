/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.Split
import StackedSVD.RankR.RMT.Stack
import StackedSVD.RMT.R1
import StackedSVD.RMT.R2
import StackedSVD.RMT.R3
import StackedSVD.RMT.R5
import StackedSVD.RMT.T

/-!
# Task U2: rank-`r` resolvent forms and the `r`-indexed `ResolventLimits`

Task U2 of `notes/archive/plan_subspacelaw.md` (sections 1.2 and 2, and the U2 row of section 3).
Everything rank one is reused verbatim from items R1, R2, R3 and T; the new content is the
pair marginals, the polarization wrappers, and the `r`-indexed interface.

1. **Marginals of the product law of U1** (plan 1.2). The split of `RankR/RMT/Split.lean`
   gives the pair `(Uᵀ Z, B)` with law `gaussianMatrix r d ⊗ gaussianMatrix p d`. The pair
   `(x_k, B)` (row `k` of the first component) has the rank-one noise law `R2.noiseLaw`
   (`measurePreserving_pairRow`), and so does the rotated pair `((x_j + s x_k)/√2, B)` for
   `s = ±1`, `j ≠ k` (`measurePreserving_pairMix`, a Householder rotation of the first
   component).
2. **Polarization wrappers** (plan 1.2 table). At complex `z` with `Im z > 0` the six forms
   per pair `(j, k)` converge: `v_jᵀ G v_k`, `g_jᵀ G g_k` to `δ_jk mᶜ(z)`, `v_jᵀ G g_k` to
   `0`, and the `G²` twins to `δ_jk (mᶜ)'(z)` and `0`. Off-diagonal terms come from the
   proved diagonal forms at the unit vectors `(v_j ± v_k)/√2` and at the rotated Gaussian
   pairs, through the complex polarization identity (`cformC_eq_polarizationC`).
3. **The interface** `ResolventLimitsR` (plan section 2): the `edge` field for `W₀` and the
   six matrix-valued limits, entrywise over `(j, k)`, at every real `z > bulkEdge c`. Item T
   transfers each complex limit to the real axis (`resolventLimitsR_of_pairLaw`). The
   derived forms `ResolventLimitsR.cform_qcol` and `.cform2_qcol` give the plan's
   `Qᵀ G₀(z) Q → (Λ + 1) m(z)` and the `G₀²` twin, entrywise.
4. **The assembly** `UnalignedModel.resolventLimitsR_of_gaussian` on the split of U1.
   The side condition that remains is `∀ N, ∑ i, n i N = r + pp N` together with
   `∀ N, 0 < pp N` (that is, `r < ∑ i, n i N` at every `N`); decision D11's tail-shift
   removal is task U7's job.

Numeric check before the proofs (a session script, `check_forms.py`, not kept): seed
`2026083002`,
`d = 400`, `M = 2`, `r = 2`, `c = 2.5`, spikes `(2.0, 1.7)`, `ψ = 0.7`, at
`x = b + 0.6` and `b + 1.6`: entrywise `|Qᵀ G₀ Q - (Λ+1)m| = 0.080, 0.051`, the `G₀²` twin
`0.064, 0.013`, every cross term at most `0.016`, split residual `2.7e-15`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

namespace FormsR

/-! ### 1. Marginals of the product law (plan 1.2) -/

section Marginals

/-- Evaluation at one row of the canonical Gaussian matrix, on plain pi types. -/
private theorem measurePreserving_rowEvalPi (rr d : ℕ) (k : Fin rr) :
    MeasurePreserving (fun A : Fin rr → Fin d → ℝ => A k)
      (Measure.pi fun _ : Fin rr => Measure.pi fun _ : Fin d => gaussianReal 0 1)
      (Measure.pi fun _ : Fin d => gaussianReal 0 1) :=
  measurePreserving_eval (μ := fun _ : Fin rr => Measure.pi fun _ : Fin d => gaussianReal 0 1) k

/-- Row `k` of a canonical Gaussian matrix is a canonical Gaussian vector. -/
theorem measurePreserving_rowEval (rr d : ℕ) (k : Fin rr) :
    MeasurePreserving (fun A : Matrix (Fin rr) (Fin d) ℝ => A k)
      (gaussianMatrix rr d) (Measure.pi fun _ : Fin d => gaussianReal 0 1) :=
  measurePreserving_rowEvalPi rr d k

/-- **The row marginal of the product law of U1** (plan 1.2): the pair `(x_k, B)` of the
rank-`r` split has the rank-one noise law of item R2. -/
theorem measurePreserving_pairRow (rr p d : ℕ) (k : Fin rr) :
    MeasurePreserving
      (fun q : Matrix (Fin rr) (Fin d) ℝ × Matrix (Fin p) (Fin d) ℝ => (q.1 k, q.2))
      ((gaussianMatrix rr d).prod (gaussianMatrix p d)) (R2.noiseLaw p d) :=
  (measurePreserving_rowEval rr d k).prod (MeasurePreserving.id _)

variable {rr : ℕ}

/-- The unit vector `(e_j + s e_k)/√2` of the polarization rotation (plan 1.2). -/
noncomputable def mixVec (rr : ℕ) (j k : Fin rr) (s : ℝ) : Fin rr → ℝ :=
  (Real.sqrt 2)⁻¹ • ((Pi.single j 1 : Fin rr → ℝ) + s • (Pi.single k 1 : Fin rr → ℝ))

private theorem sqrt_two_inv_sq : (Real.sqrt 2)⁻¹ * (Real.sqrt 2)⁻¹ = (2 : ℝ)⁻¹ := by
  rw [← mul_inv, Real.mul_self_sqrt (by norm_num : (0 : ℝ) ≤ 2)]

theorem mixVec_dotProduct_self {j k : Fin rr} (hjk : j ≠ k) {s : ℝ}
    (hs : s = 1 ∨ s = -1) : mixVec rr j k s ⬝ᵥ mixVec rr j k s = 1 := by
  have hjj : (Pi.single j 1 : Fin rr → ℝ) ⬝ᵥ Pi.single j 1 = 1 := by
    rw [single_dotProduct, Pi.single_eq_same, mul_one]
  have hkk : (Pi.single k 1 : Fin rr → ℝ) ⬝ᵥ Pi.single k 1 = 1 := by
    rw [single_dotProduct, Pi.single_eq_same, mul_one]
  have hjk0 : (Pi.single j 1 : Fin rr → ℝ) ⬝ᵥ Pi.single k 1 = 0 := by
    rw [single_dotProduct, one_mul, Pi.single_eq_of_ne hjk]
  have hkj0 : (Pi.single k 1 : Fin rr → ℝ) ⬝ᵥ Pi.single j 1 = 0 := by
    rw [single_dotProduct, one_mul, Pi.single_eq_of_ne (Ne.symm hjk)]
  have hss : s * s = 1 := by rcases hs with h | h <;> rw [h] <;> norm_num
  rw [mixVec]
  have hexp : ((Real.sqrt 2)⁻¹ • ((Pi.single j 1 : Fin rr → ℝ) + s • Pi.single k 1)) ⬝ᵥ
      ((Real.sqrt 2)⁻¹ • ((Pi.single j 1 : Fin rr → ℝ) + s • Pi.single k 1))
      = (Real.sqrt 2)⁻¹ * (Real.sqrt 2)⁻¹ *
        (((Pi.single j 1 : Fin rr → ℝ) ⬝ᵥ Pi.single j 1)
          + s * ((Pi.single j 1 : Fin rr → ℝ) ⬝ᵥ Pi.single k 1)
          + s * ((Pi.single k 1 : Fin rr → ℝ) ⬝ᵥ Pi.single j 1)
          + s * s * ((Pi.single k 1 : Fin rr → ℝ) ⬝ᵥ Pi.single k 1)) := by
    simp only [smul_dotProduct, dotProduct_smul, add_dotProduct, dotProduct_add, smul_eq_mul]
    ring
  rw [hexp, hjj, hkk, hjk0, hkj0, sqrt_two_inv_sq, hss]
  norm_num

theorem single_ne_mixVec {j k : Fin rr} (hjk : j ≠ k) (s : ℝ) :
    (Pi.single j 1 : Fin rr → ℝ) ≠ mixVec rr j k s := by
  intro h
  have hj := congrFun h j
  rw [Pi.single_eq_same] at hj
  have hmix : mixVec rr j k s j = (Real.sqrt 2)⁻¹ := by
    simp only [mixVec, Pi.smul_apply, Pi.add_apply, smul_eq_mul, Pi.single_eq_same,
      Pi.single_eq_of_ne hjk, mul_zero, add_zero, mul_one]
  have h12 : (1 : ℝ) < Real.sqrt 2 := by
    have : Real.sqrt 1 < Real.sqrt 2 := Real.sqrt_lt_sqrt (by norm_num) (by norm_num)
    simpa using this
  have hs2 : Real.sqrt 2 ≠ 0 := by positivity
  rw [hmix] at hj
  have hone : Real.sqrt 2 = 1 := by
    have h2 : Real.sqrt 2 * (Real.sqrt 2)⁻¹ = 1 := mul_inv_cancel₀ hs2
    rw [← hj, mul_one] at h2
    exact h2
  linarith

/-- Row `j` of `householder (e_j - mixVec) * A` is the mixed row `(A_j + s A_k)/√2`. -/
private theorem householder_mix_row {d : ℕ} {j k : Fin rr} (hjk : j ≠ k) {s : ℝ}
    (hs : s = 1 ∨ s = -1) (A : Matrix (Fin rr) (Fin d) ℝ) (l : Fin d) :
    (householder ((Pi.single j 1 : Fin rr → ℝ) - mixVec rr j k s) * A) j l
      = (Real.sqrt 2)⁻¹ * (A j l + s * A k l) := by
  set w : Fin rr → ℝ := (Pi.single j 1 : Fin rr → ℝ) - mixVec rr j k s with hwdef
  set H : Matrix (Fin rr) (Fin rr) ℝ := householder w with hHdef
  have hone : (Pi.single j 1 : Fin rr → ℝ) ⬝ᵥ Pi.single j 1 = 1 := by
    rw [single_dotProduct, Pi.single_eq_same, mul_one]
  have hcol : H *ᵥ (Pi.single j 1 : Fin rr → ℝ) = mixVec rr j k s :=
    householder_sub_apply hone (mixVec_dotProduct_self hjk hs) (single_ne_mixVec hjk s)
  have hsym : Hᵀ = H := householder_transpose w
  have hHji : ∀ i : Fin rr, H j i = mixVec rr j k s i := by
    intro i
    have h1 : H j i = H i j := congrFun (congrFun hsym i) j
    rw [h1, ← hcol]
    have h3 := dotProduct_single (v := fun l => H i l) (1 : ℝ) j
    simpa [Matrix.mulVec, dotProduct] using h3.symm
  have hsumj : ∑ i, (Pi.single j 1 : Fin rr → ℝ) i * A i l = A j l := by
    have h3 := single_dotProduct (v := fun i => A i l) (1 : ℝ) j
    simpa [dotProduct] using h3
  have hsumk : ∑ i, (Pi.single k 1 : Fin rr → ℝ) i * A i l = A k l := by
    have h3 := single_dotProduct (v := fun i => A i l) (1 : ℝ) k
    simpa [dotProduct] using h3
  have hexp : ∀ i : Fin rr, mixVec rr j k s i * A i l
      = (Real.sqrt 2)⁻¹ * ((Pi.single j 1 : Fin rr → ℝ) i * A i l)
        + (Real.sqrt 2)⁻¹ * s * ((Pi.single k 1 : Fin rr → ℝ) i * A i l) := by
    intro i
    simp only [mixVec, Pi.smul_apply, Pi.add_apply, smul_eq_mul]
    ring
  calc (H * A) j l = ∑ i, H j i * A i l := Matrix.mul_apply
    _ = ∑ i, mixVec rr j k s i * A i l :=
        Finset.sum_congr rfl fun i _ => by rw [hHji i]
    _ = (Real.sqrt 2)⁻¹ * (∑ i, (Pi.single j 1 : Fin rr → ℝ) i * A i l)
        + (Real.sqrt 2)⁻¹ * s * (∑ i, (Pi.single k 1 : Fin rr → ℝ) i * A i l) := by
        rw [Finset.sum_congr rfl fun i _ => hexp i, Finset.sum_add_distrib,
          ← Finset.mul_sum, ← Finset.mul_sum]
    _ = (Real.sqrt 2)⁻¹ * (A j l + s * A k l) := by rw [hsumj, hsumk]; ring

/-- **The rotated marginal of the product law of U1** (plan 1.2): the pair
`((x_j + s x_k)/√2, B)` with `s = ±1` and `j ≠ k` has the rank-one noise law. The rotation
is the Householder reflection that sends `e_j` to `(e_j + s e_k)/√2`, a left rotation of
the first component (`measurePreserving_mul_left`). -/
theorem measurePreserving_pairMix (rr p d : ℕ) {j k : Fin rr} (hjk : j ≠ k) {s : ℝ}
    (hs : s = 1 ∨ s = -1) :
    MeasurePreserving
      (fun q : Matrix (Fin rr) (Fin d) ℝ × Matrix (Fin p) (Fin d) ℝ =>
        ((fun l => (Real.sqrt 2)⁻¹ * (q.1 j l + s * q.1 k l)), q.2))
      ((gaussianMatrix rr d).prod (gaussianMatrix p d)) (R2.noiseLaw p d) := by
  set w : Fin rr → ℝ := (Pi.single j 1 : Fin rr → ℝ) - mixVec rr j k s with hwdef
  have hne : (Pi.single j 1 : Fin rr → ℝ) ≠ mixVec rr j k s := single_ne_mixVec hjk s
  have hw0 : w ≠ 0 := sub_ne_zero.mpr hne
  have hww : w ⬝ᵥ w ≠ 0 := fun h => hw0 (dotProduct_self_eq_zero.mp h)
  have hH : (householder w)ᵀ * householder w = 1 := householder_orth hww
  have hcomp := (measurePreserving_pairRow rr p d j).comp
    ((measurePreserving_mul_left hH d).prod
      (MeasurePreserving.id (gaussianMatrix p d)))
  have hfun : ((fun q : Matrix (Fin rr) (Fin d) ℝ × Matrix (Fin p) (Fin d) ℝ =>
        (q.1 j, q.2)) ∘
        Prod.map (fun A : Matrix (Fin rr) (Fin d) ℝ => householder w * A) id)
      = fun q : Matrix (Fin rr) (Fin d) ℝ × Matrix (Fin p) (Fin d) ℝ =>
        ((fun l => (Real.sqrt 2)⁻¹ * (q.1 j l + s * q.1 k l)), q.2) := by
    funext q
    have hrow : (householder w * q.1) j
        = fun l => (Real.sqrt 2)⁻¹ * (q.1 j l + s * q.1 k l) :=
      funext fun l => householder_mix_row hjk hs q.1 l
    change ((householder w * q.1) j, q.2) = _
    rw [hrow]
  rwa [hfun] at hcomp

end Marginals

/-! ### 2. Complex bilinearity, symmetry, and polarization of the forms -/

section ComplexForms

variable {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ} {z : ℂ}

theorem cvec_add (x y : Fin D → ℝ) : R4C.cvec (x + y) = R4C.cvec x + R4C.cvec y := by
  funext a
  simp [R4C.cvec]

theorem cformC_add_left (W : Matrix (Fin D) (Fin D) ℝ) (z : ℂ) (x x' y : Fin D → ℝ) :
    R4C.cformC W z (x + x') y = R4C.cformC W z x y + R4C.cformC W z x' y := by
  simp only [R4C.cformC, cvec_add, add_dotProduct]

theorem cformC_add_right (W : Matrix (Fin D) (Fin D) ℝ) (z : ℂ) (x y y' : Fin D → ℝ) :
    R4C.cformC W z x (y + y') = R4C.cformC W z x y + R4C.cformC W z x y' := by
  simp only [R4C.cformC, cvec_add, Matrix.mulVec_add, dotProduct_add]

theorem cform2C_add_left (W : Matrix (Fin D) (Fin D) ℝ) (z : ℂ) (x x' y : Fin D → ℝ) :
    R4C.cform2C W z (x + x') y = R4C.cform2C W z x y + R4C.cform2C W z x' y := by
  simp only [R4C.cform2C, cvec_add, add_dotProduct]

theorem cform2C_add_right (W : Matrix (Fin D) (Fin D) ℝ) (z : ℂ) (x y y' : Fin D → ℝ) :
    R4C.cform2C W z x (y + y') = R4C.cform2C W z x y + R4C.cform2C W z x y' := by
  simp only [R4C.cform2C, cvec_add, Matrix.mulVec_add, dotProduct_add]

/-- The complex resolvent form is symmetric (the kernel is complex symmetric). -/
theorem cformC_comm (hW : W.IsHermitian) (hz : z.im ≠ 0) (x y : Fin D → ℝ) :
    R4C.cformC W z x y = R4C.cformC W z y x := by
  rw [R4C.cformC_eq_sum hW hz, R4C.cformC_eq_sum hW hz]
  exact Finset.sum_congr rfl fun a _ => by ring_nf

theorem cform2C_comm (hW : W.IsHermitian) (hz : z.im ≠ 0) (x y : Fin D → ℝ) :
    R4C.cform2C W z x y = R4C.cform2C W z y x := by
  rw [R4C.cform2C_eq_sum hW hz, R4C.cform2C_eq_sum hW hz]
  exact Finset.sum_congr rfl fun a _ => by ring_nf

/-- **Complex polarization**: the cross form from three quadratic forms. This is the rank-`r`
polarization device of plan 1.2, at complex `z`. -/
theorem cformC_eq_polarizationC (hW : W.IsHermitian) (hz : z.im ≠ 0) (x y : Fin D → ℝ) :
    R4C.cformC W z x y
      = (R4C.qformC W z (x + y) - R4C.qformC W z x - R4C.qformC W z y) / 2 := by
  have hexp : R4C.qformC W z (x + y)
      = R4C.qformC W z x + R4C.cformC W z x y + R4C.cformC W z y x + R4C.qformC W z y := by
    simp only [R4C.qformC, cformC_add_left, cformC_add_right]
    ring
  rw [hexp, cformC_comm hW hz y x]
  ring

theorem cform2C_eq_polarizationC (hW : W.IsHermitian) (hz : z.im ≠ 0) (x y : Fin D → ℝ) :
    R4C.cform2C W z x y
      = (R4C.qform2C W z (x + y) - R4C.qform2C W z x - R4C.qform2C W z y) / 2 := by
  have hexp : R4C.qform2C W z (x + y)
      = R4C.qform2C W z x + R4C.cform2C W z x y + R4C.cform2C W z y x
        + R4C.qform2C W z y := by
    simp only [R4C.qform2C, cform2C_add_left, cform2C_add_right]
    ring
  rw [hexp, cform2C_comm hW hz y x]
  ring

/-- `gOf` of the mixed row is the mix of the `gOf`s. -/
theorem gOf_mixRow {D : ℕ} (x y : Fin D → ℝ) (s : ℝ) :
    R2.gOf (fun l => (Real.sqrt 2)⁻¹ * (x l + s * y l))
      = (Real.sqrt 2)⁻¹ • (R2.gOf x + s • R2.gOf y) := by
  funext a
  simp only [R2.gOf, Pi.smul_apply, Pi.add_apply, smul_eq_mul]
  ring

/-- Scaling the mix back: `g_j + s g_k = √2 • gOf(mixed row)`. -/
theorem add_smul_gOf_eq (x y : Fin D → ℝ) (s : ℝ) :
    R2.gOf x + s • R2.gOf y
      = Real.sqrt 2 • R2.gOf (fun l => (Real.sqrt 2)⁻¹ * (x l + s * y l)) := by
  rw [gOf_mixRow, smul_smul, mul_inv_cancel₀ (Real.sqrt_ne_zero'.mpr (by norm_num)),
    one_smul]

theorem qformC_eq_cformC {D : ℕ} (W : Matrix (Fin D) (Fin D) ℝ) (z : ℂ) (y : Fin D → ℝ) :
    R4C.qformC W z y = R4C.cformC W z y y := rfl

theorem qform2C_eq_cform2C {D : ℕ} (W : Matrix (Fin D) (Fin D) ℝ) (z : ℂ) (y : Fin D → ℝ) :
    R4C.qform2C W z y = R4C.cform2C W z y y := rfl

end ComplexForms

/-! ### 3. The complex limits per pair (plan 1.2 table) -/

section ComplexLimits

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {pN dN : ℕ → ℕ} {r : ℕ} {c : ℝ} {z : ℂ}
  (ZZ : ∀ N, Ω N → Matrix (Fin r) (Fin (dN N)) ℝ × Matrix (Fin (pN N)) (Fin (dN N)) ℝ)

/-- The law of the pair `(x_k, B)` read off the joint pair of U1. -/
theorem hasLaw_pairRow
    (hZZ : ∀ N, HasLaw (ZZ N)
      ((gaussianMatrix r (dN N)).prod (gaussianMatrix (pN N) (dN N))) (μ N)) (k : Fin r) :
    ∀ N, HasLaw (fun ω => ((ZZ N ω).1 k, (ZZ N ω).2)) (R2.noiseLaw (pN N) (dN N)) (μ N) :=
  fun N => (measurePreserving_pairRow r (pN N) (dN N) k).fun_comp_hasLaw (hZZ N)

/-- The law of the rotated pair `((x_j + s x_k)/√2, B)` read off the joint pair of U1. -/
theorem hasLaw_pairMix
    (hZZ : ∀ N, HasLaw (ZZ N)
      ((gaussianMatrix r (dN N)).prod (gaussianMatrix (pN N) (dN N))) (μ N))
    {j k : Fin r} (hjk : j ≠ k) {s : ℝ} (hs : s = 1 ∨ s = -1) :
    ∀ N, HasLaw (fun ω =>
      ((fun l => (Real.sqrt 2)⁻¹ * ((ZZ N ω).1 j l + s * (ZZ N ω).1 k l)), (ZZ N ω).2))
      (R2.noiseLaw (pN N) (dN N)) (μ N) :=
  fun N => (measurePreserving_pairMix r (pN N) (dN N) hjk hs).fun_comp_hasLaw (hZZ N)

private theorem norm_two_C : ‖(2 : ℂ)‖ = 2 := by
  have h : ((2 : ℝ) : ℂ) = (2 : ℂ) := by norm_num
  rw [← h, Complex.norm_real, Real.norm_eq_abs]
  norm_num

private theorem norm_eq_one_of_inner {D : ℕ} {x : EuclideanSpace ℝ (Fin D)}
    (h : ⟪x, x⟫_ℝ = 1) : ‖x‖ = 1 := by
  have hsq : ‖x‖ ^ 2 = 1 := by rw [← real_inner_self_eq_norm_sq]; exact h
  have hx : ‖x‖ = Real.sqrt (‖x‖ ^ 2) := (Real.sqrt_sq (norm_nonneg x)).symm
  rw [hx, hsq, Real.sqrt_one]

/-- **The polarization step at complex `z`** (plan 1.2): diagonal limits at `a`, `b` and at
the mix `u` with `a + b = √2 u` give the cross limit `0`. First order. -/
theorem tendstoInProb_cformC_cross
    {Wf : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ}
    (hWf : ∀ N ω, (Wf N ω).IsHermitian) (hz : 0 < z.im)
    {af bf uf : (N : ℕ) → Ω N → Fin (dN N) → ℝ} {L : ℂ}
    (hsum : ∀ N ω, af N ω + bf N ω = Real.sqrt 2 • uf N ω)
    (hqa : TendstoInProb μ (fun N ω => ‖R4C.qformC (Wf N ω) z (af N ω) - L‖) 0)
    (hqb : TendstoInProb μ (fun N ω => ‖R4C.qformC (Wf N ω) z (bf N ω) - L‖) 0)
    (hqu : TendstoInProb μ (fun N ω => ‖R4C.qformC (Wf N ω) z (uf N ω) - L‖) 0) :
    TendstoInProb μ (fun N ω => ‖R4C.cformC (Wf N ω) z (af N ω) (bf N ω)‖) 0 := by
  have h2C : ((Real.sqrt 2 ^ 2 : ℝ) : ℂ) = 2 := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 2)]
    norm_num
  refine TendstoInProb.of_le
    (g := fun N ω => ‖R4C.qformC (Wf N ω) z (uf N ω) - L‖
      + (‖R4C.qformC (Wf N ω) z (af N ω) - L‖ + ‖R4C.qformC (Wf N ω) z (bf N ω) - L‖) * 2⁻¹)
    (fun N => ?_) ?_
  · filter_upwards with ω
    rw [sub_zero, abs_of_nonneg (norm_nonneg _)]
    have hpol := cformC_eq_polarizationC (hWf N ω) hz.ne' (af N ω) (bf N ω)
    have hkey : R4C.qformC (Wf N ω) z (af N ω + bf N ω)
        = 2 * R4C.qformC (Wf N ω) z (uf N ω) := by
      rw [hsum N ω, R2.qformC_smul, h2C]
    have heq : R4C.cformC (Wf N ω) z (af N ω) (bf N ω)
        = ((2 : ℂ) * (R4C.qformC (Wf N ω) z (uf N ω) - L)
          - (R4C.qformC (Wf N ω) z (af N ω) - L)
          - (R4C.qformC (Wf N ω) z (bf N ω) - L)) / 2 := by
      rw [hpol, hkey]
      ring
    rw [heq, norm_div, norm_two_C]
    have h4 : ‖(2 : ℂ) * (R4C.qformC (Wf N ω) z (uf N ω) - L)‖
        = 2 * ‖R4C.qformC (Wf N ω) z (uf N ω) - L‖ := by
      rw [norm_mul, norm_two_C]
    have h5 := norm_sub_le ((2 : ℂ) * (R4C.qformC (Wf N ω) z (uf N ω) - L)
        - (R4C.qformC (Wf N ω) z (af N ω) - L)) (R4C.qformC (Wf N ω) z (bf N ω) - L)
    have h6 := norm_sub_le ((2 : ℂ) * (R4C.qformC (Wf N ω) z (uf N ω) - L))
      (R4C.qformC (Wf N ω) z (af N ω) - L)
    rw [h4] at h6
    linarith
  · have hlim := hqu.add ((hqa.add hqb).mul_const 2⁻¹)
    simpa using hlim

/-- The polarization step, second order. -/
theorem tendstoInProb_cform2C_cross
    {Wf : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ}
    (hWf : ∀ N ω, (Wf N ω).IsHermitian) (hz : 0 < z.im)
    {af bf uf : (N : ℕ) → Ω N → Fin (dN N) → ℝ} {L : ℂ}
    (hsum : ∀ N ω, af N ω + bf N ω = Real.sqrt 2 • uf N ω)
    (hqa : TendstoInProb μ (fun N ω => ‖R4C.qform2C (Wf N ω) z (af N ω) - L‖) 0)
    (hqb : TendstoInProb μ (fun N ω => ‖R4C.qform2C (Wf N ω) z (bf N ω) - L‖) 0)
    (hqu : TendstoInProb μ (fun N ω => ‖R4C.qform2C (Wf N ω) z (uf N ω) - L‖) 0) :
    TendstoInProb μ (fun N ω => ‖R4C.cform2C (Wf N ω) z (af N ω) (bf N ω)‖) 0 := by
  have h2C : ((Real.sqrt 2 ^ 2 : ℝ) : ℂ) = 2 := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 2)]
    norm_num
  refine TendstoInProb.of_le
    (g := fun N ω => ‖R4C.qform2C (Wf N ω) z (uf N ω) - L‖
      + (‖R4C.qform2C (Wf N ω) z (af N ω) - L‖
        + ‖R4C.qform2C (Wf N ω) z (bf N ω) - L‖) * 2⁻¹)
    (fun N => ?_) ?_
  · filter_upwards with ω
    rw [sub_zero, abs_of_nonneg (norm_nonneg _)]
    have hpol := cform2C_eq_polarizationC (hWf N ω) hz.ne' (af N ω) (bf N ω)
    have hkey : R4C.qform2C (Wf N ω) z (af N ω + bf N ω)
        = 2 * R4C.qform2C (Wf N ω) z (uf N ω) := by
      rw [hsum N ω, R2.qform2C_smul, h2C]
    have heq : R4C.cform2C (Wf N ω) z (af N ω) (bf N ω)
        = ((2 : ℂ) * (R4C.qform2C (Wf N ω) z (uf N ω) - L)
          - (R4C.qform2C (Wf N ω) z (af N ω) - L)
          - (R4C.qform2C (Wf N ω) z (bf N ω) - L)) / 2 := by
      rw [hpol, hkey]
      ring
    rw [heq, norm_div, norm_two_C]
    have h4 : ‖(2 : ℂ) * (R4C.qform2C (Wf N ω) z (uf N ω) - L)‖
        = 2 * ‖R4C.qform2C (Wf N ω) z (uf N ω) - L‖ := by
      rw [norm_mul, norm_two_C]
    have h5 := norm_sub_le ((2 : ℂ) * (R4C.qform2C (Wf N ω) z (uf N ω) - L)
        - (R4C.qform2C (Wf N ω) z (af N ω) - L)) (R4C.qform2C (Wf N ω) z (bf N ω) - L)
    have h6 := norm_sub_le ((2 : ℂ) * (R4C.qform2C (Wf N ω) z (uf N ω) - L))
      (R4C.qform2C (Wf N ω) z (af N ω) - L)
    rw [h4] at h6
    linarith
  · have hlim := hqu.add ((hqa.add hqb).mul_const 2⁻¹)
    simpa using hlim

/-! #### The diagonal `g`-forms through the row marginal -/

/-- `g_kᵀ G(z) g_k → mᶜ(z)`: item R2 at the pair `(x_k, B)`, centered by R1a. -/
theorem tendstoInProb_qformC_g (hz : 0 < z.im)
    (hZZ : ∀ N, HasLaw (ZZ N)
      ((gaussianMatrix r (dN N)).prod (gaussianMatrix (pN N) (dN N))) (μ N))
    (hd : Tendsto dN atTop atTop) (k : Fin r)
    (hR1a : TendstoInProb μ
      (fun N ω => ‖R4C.stieltjesC (R2.W0 (ZZ N ω).2) z - MP.mC c z‖) 0) :
    TendstoInProb μ (fun N ω =>
      ‖R4C.qformC (R2.W0 (ZZ N ω).2) z (R2.gOf ((ZZ N ω).1 k)) - MP.mC c z‖) 0 :=
  R2.tendstoInProb_norm_sub_trans
    (R2.tendstoInProb_qformC_gOf_sub (fun N ω => ((ZZ N ω).1 k, (ZZ N ω).2)) hz
      (hasLaw_pairRow ZZ hZZ k) hd) hR1a

/-- `g_kᵀ G(z)² g_k → (mᶜ)'(z)`. -/
theorem tendstoInProb_qform2C_g (hz : 0 < z.im)
    (hZZ : ∀ N, HasLaw (ZZ N)
      ((gaussianMatrix r (dN N)).prod (gaussianMatrix (pN N) (dN N))) (μ N))
    (hd : Tendsto dN atTop atTop) (k : Fin r)
    (hR1b : TendstoInProb μ
      (fun N ω => ‖R4C.stieltjes2C (R2.W0 (ZZ N ω).2) z - MP.mCDeriv c z‖) 0) :
    TendstoInProb μ (fun N ω =>
      ‖R4C.qform2C (R2.W0 (ZZ N ω).2) z (R2.gOf ((ZZ N ω).1 k)) - MP.mCDeriv c z‖) 0 :=
  R2.tendstoInProb_norm_sub_trans
    (R2.tendstoInProb_qform2C_gOf_sub (fun N ω => ((ZZ N ω).1 k, (ZZ N ω).2)) hz
      (hasLaw_pairRow ZZ hZZ k) hd) hR1b

/-- The `g`-form at the rotated pair, first order (the `u` input of polarization). -/
theorem tendstoInProb_qformC_gmix (hz : 0 < z.im)
    (hZZ : ∀ N, HasLaw (ZZ N)
      ((gaussianMatrix r (dN N)).prod (gaussianMatrix (pN N) (dN N))) (μ N))
    (hd : Tendsto dN atTop atTop) {j k : Fin r} (hjk : j ≠ k)
    (hR1a : TendstoInProb μ
      (fun N ω => ‖R4C.stieltjesC (R2.W0 (ZZ N ω).2) z - MP.mC c z‖) 0) :
    TendstoInProb μ (fun N ω => ‖R4C.qformC (R2.W0 (ZZ N ω).2) z
      (R2.gOf (fun l => (Real.sqrt 2)⁻¹ * ((ZZ N ω).1 j l + 1 * (ZZ N ω).1 k l)))
      - MP.mC c z‖) 0 :=
  R2.tendstoInProb_norm_sub_trans
    (R2.tendstoInProb_qformC_gOf_sub (fun N ω =>
        ((fun l => (Real.sqrt 2)⁻¹ * ((ZZ N ω).1 j l + 1 * (ZZ N ω).1 k l)), (ZZ N ω).2))
      hz (hasLaw_pairMix ZZ hZZ hjk (Or.inl rfl)) hd) hR1a

/-- The `g`-form at the rotated pair, second order. -/
theorem tendstoInProb_qform2C_gmix (hz : 0 < z.im)
    (hZZ : ∀ N, HasLaw (ZZ N)
      ((gaussianMatrix r (dN N)).prod (gaussianMatrix (pN N) (dN N))) (μ N))
    (hd : Tendsto dN atTop atTop) {j k : Fin r} (hjk : j ≠ k)
    (hR1b : TendstoInProb μ
      (fun N ω => ‖R4C.stieltjes2C (R2.W0 (ZZ N ω).2) z - MP.mCDeriv c z‖) 0) :
    TendstoInProb μ (fun N ω => ‖R4C.qform2C (R2.W0 (ZZ N ω).2) z
      (R2.gOf (fun l => (Real.sqrt 2)⁻¹ * ((ZZ N ω).1 j l + 1 * (ZZ N ω).1 k l)))
      - MP.mCDeriv c z‖) 0 :=
  R2.tendstoInProb_norm_sub_trans
    (R2.tendstoInProb_qform2C_gOf_sub (fun N ω =>
        ((fun l => (Real.sqrt 2)⁻¹ * ((ZZ N ω).1 j l + 1 * (ZZ N ω).1 k l)), (ZZ N ω).2))
      hz (hasLaw_pairMix ZZ hZZ hjk (Or.inl rfl)) hd) hR1b

/-- The rank-`r` split sum identity: `g_j + g_k = √2 gOf(mixed row)`. -/
theorem gOf_add_eq_sqrt_two_smul {D : ℕ} (x y : Fin D → ℝ) :
    R2.gOf x + R2.gOf y
      = Real.sqrt 2 • R2.gOf (fun l => (Real.sqrt 2)⁻¹ * (x l + 1 * y l)) := by
  have h := add_smul_gOf_eq x y 1
  rwa [one_smul] at h

end ComplexLimits

end FormsR

/-! ### 4. The `r`-indexed interface (plan section 2, U2 row) -/

section Interface

open R4

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]

/-- **(H1) and (H2) at rank `r`** (plan 1.2 and section 2): the `edge` bound for the block
`W₀` and the six resolvent limits, entrywise over the pair `(j, k)`, at every real
`z > bulkEdge c`. This is the rank-`r` twin of `SpikedModel.ResolventLimits` (`RMT/R5.lean`,
decision D6). The deterministic directions are the `r` spike vectors `v j N` (orthonormal for
each `N`); the random directions are the scaled Gaussian columns `g j N ω` of the split of
`RankR/RMT/Split.lean`. Entry `(j, k)` of the plan's matrix limit `Qᵀ G₀(z) Q → (Λ + 1) m(z)`
follows by bilinearity (`ResolventLimitsR.cform_qcol`). -/
structure ResolventLimitsR (μ : ∀ N, Measure (Ω N)) {dN : ℕ → ℕ} {r : ℕ}
    (W₀ : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ)
    (hsymm : ∀ N ω, (W₀ N ω).IsHermitian)
    (v : Fin r → (N : ℕ) → EuclideanSpace ℝ (Fin (dN N)))
    (g : Fin r → (N : ℕ) → Ω N → Fin (dN N) → ℝ) (c : ℝ) : Prop where
  /-- (H1), from item R3 (`tendsto_measure_lamMax_le`) on the block of the rank-`r` split. -/
  edge : ∀ ε > 0, Tendsto (fun N => μ N
    {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ bulkEdge c + ε}) atTop (𝓝 1)
  /-- (H2), `v_jᵀ G₀ v_k → δ_jk m`. -/
  vv : ∀ j k : Fin r, ∀ z, bulkEdge c < z → TendstoInProb μ
    (fun N ω => cform (W₀ N ω) z (WithLp.ofLp (v j N)) (WithLp.ofLp (v k N)))
    (if j = k then MP.m c z else 0)
  /-- (H2), `g_jᵀ G₀ g_k → δ_jk m`. -/
  gg : ∀ j k : Fin r, ∀ z, bulkEdge c < z → TendstoInProb μ
    (fun N ω => cform (W₀ N ω) z (g j N ω) (g k N ω)) (if j = k then MP.m c z else 0)
  /-- (H2), `v_jᵀ G₀ g_k → 0`, every pair. -/
  vg : ∀ j k : Fin r, ∀ z, bulkEdge c < z → TendstoInProb μ
    (fun N ω => cform (W₀ N ω) z (WithLp.ofLp (v j N)) (g k N ω)) 0
  /-- (H2), `v_jᵀ G₀² v_k → δ_jk m'`. -/
  vv2 : ∀ j k : Fin r, ∀ z, bulkEdge c < z → TendstoInProb μ
    (fun N ω => cform2 (W₀ N ω) z (WithLp.ofLp (v j N)) (WithLp.ofLp (v k N)))
    (if j = k then MP.mDeriv c z else 0)
  /-- (H2), `g_jᵀ G₀² g_k → δ_jk m'`. -/
  gg2 : ∀ j k : Fin r, ∀ z, bulkEdge c < z → TendstoInProb μ
    (fun N ω => cform2 (W₀ N ω) z (g j N ω) (g k N ω)) (if j = k then MP.mDeriv c z else 0)
  /-- (H2), `v_jᵀ G₀² g_k → 0`, every pair. -/
  vg2 : ∀ j k : Fin r, ∀ z, bulkEdge c < z → TendstoInProb μ
    (fun N ω => cform2 (W₀ N ω) z (WithLp.ofLp (v j N)) (g k N ω)) 0

end Interface

/-! ### 5. The transfer to the real axis: the interface from the pair law -/

section Transfer

open R4

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {pN dN : ℕ → ℕ} {r : ℕ}

/-- **The `r`-indexed interface from the pair law** (plan 1.2, the T wrappers). Model-free:
the data are a block family `W₀` (pointwise `R2.W0` of the second pair component), the
orthonormal spike families `v`, the scaled Gaussian columns `g` (pointwise `R2.gOf` of the
rows of the first pair component), and the edge bound in both shapes. Items R1 and R2 give
the complex limits per pair through the marginals of section 3, and item T transfers each
one to every real `z > bulkEdge c`. -/
theorem resolventLimitsR_of_pairLaw [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (hd : Tendsto dN atTop atTop) (hp : ∀ N, 0 < pN N)
    (hcN : Tendsto (fun N => (pN N : ℝ) / dN N) atTop (𝓝 c))
    (W₀ : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ)
    (hsymm : ∀ N ω, (W₀ N ω).IsHermitian)
    (v : Fin r → (N : ℕ) → EuclideanSpace ℝ (Fin (dN N)))
    (g : Fin r → (N : ℕ) → Ω N → Fin (dN N) → ℝ)
    (ZZ : ∀ N, Ω N → Matrix (Fin r) (Fin (dN N)) ℝ × Matrix (Fin (pN N)) (Fin (dN N)) ℝ)
    (hZZ : ∀ N, HasLaw (ZZ N)
      ((gaussianMatrix r (dN N)).prod (gaussianMatrix (pN N) (dN N))) (μ N))
    (hW0eq : ∀ N ω, R2.W0 (ZZ N ω).2 = W₀ N ω)
    (hgeq : ∀ (k : Fin r) (N : ℕ) (ω : Ω N), R2.gOf ((ZZ N ω).1 k) = g k N ω)
    (hv : ∀ (a b : Fin r) (N : ℕ), ⟪v a N, v b N⟫_ℝ = if a = b then 1 else 0)
    (hedge : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ bulkEdge c + ε}) atTop (𝓝 1))
    (hedgeC : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (W₀ N ω) (hsymm N ω) ≤ bulkEdge c + ε}ᶜ) atTop (𝓝 0)) :
    ResolventLimitsR μ W₀ hsymm v g c := by
  classical
  have hBlaw : ∀ N, HasLaw (fun ω => (ZZ N ω).2) (gaussianMatrix (pN N) (dN N)) (μ N) :=
    fun N => measurePreserving_snd.fun_comp_hasLaw (hZZ N)
  have hR1a : ∀ z : ℂ, 0 < z.im → TendstoInProb μ
      (fun N ω => ‖R4C.stieltjesC (R2.W0 (ZZ N ω).2) z - MP.mC c z‖) 0 := fun z hz =>
    R1.tendstoInProb_stieltjesC hc hz hd hp hcN (fun N ω => (ZZ N ω).2) hBlaw
  have hR1b : ∀ z : ℂ, 0 < z.im → TendstoInProb μ
      (fun N ω => ‖R4C.stieltjes2C (R2.W0 (ZZ N ω).2) z - MP.mCDeriv c z‖) 0 := fun z hz =>
    R1.tendstoInProb_stieltjes2C hc hz hd hp hcN (fun N ω => (ZZ N ω).2) hBlaw
  have hnorm : ∀ (a : Fin r) (N : ℕ), ‖v a N‖ = 1 := fun a N =>
    FormsR.norm_eq_one_of_inner (by have h := hv a a N; rwa [if_pos rfl] at h)
  have hvdot : ∀ (a : Fin r) (N : ℕ),
      WithLp.ofLp (v a N) ⬝ᵥ WithLp.ofLp (v a N) = 1 :=
    fun a N => R2.dotProduct_ofLp_self (hnorm a N)
  have hsqrt2 : Real.sqrt 2 ≠ 0 := Real.sqrt_ne_zero'.mpr (by norm_num)
  -- the unit mix family of the deterministic directions
  have hup : ∀ j k : Fin r, j ≠ k → ∀ N : ℕ,
      ‖(Real.sqrt 2)⁻¹ • (v j N + v k N)‖ = 1 := by
    intro j k hjk N
    have hin : ⟪v j N, v k N⟫_ℝ = 0 := by
      have h := hv j k N
      rwa [if_neg hjk] at h
    have hsq : ‖v j N + v k N‖ ^ 2 = 2 := by
      rw [norm_add_sq_real, hnorm j N, hnorm k N, hin]
      norm_num
    have hnn : ‖v j N + v k N‖ = Real.sqrt 2 := by
      rw [← Real.sqrt_sq (norm_nonneg (v j N + v k N)), hsq]
    rw [norm_smul, hnn, norm_inv, Real.norm_eq_abs, abs_of_nonneg (Real.sqrt_nonneg 2),
      inv_mul_cancel₀ hsqrt2]
  have hvsum : ∀ (j k : Fin r) (N : ℕ),
      WithLp.ofLp (v j N) + WithLp.ofLp (v k N)
        = Real.sqrt 2 • WithLp.ofLp ((Real.sqrt 2)⁻¹ • (v j N + v k N)) := by
    intro j k N
    rw [WithLp.ofLp_smul, WithLp.ofLp_add, smul_smul, mul_inv_cancel₀ hsqrt2, one_smul]
  -- the norm events of item T
  have hgdot : ∀ k : Fin r, TendstoInProb μ (fun N ω => g k N ω ⬝ᵥ g k N ω) 1 := by
    intro k
    have h := R2.tendstoInProb_dotProduct_gOf (fun N ω => ((ZZ N ω).1 k, (ZZ N ω).2))
      (FormsR.hasLaw_pairRow ZZ hZZ k) hd
    simpa only [hgeq] using h
  have hgbad : ∀ k : Fin r,
      Tendsto (fun N => μ N {ω | g k N ω ⬝ᵥ g k N ω ≤ 4}ᶜ) atTop (𝓝 0) := by
    intro k
    refine tendsto_measure_zero_of_subset
      (t := fun N => {ω | (3 : ℝ) ≤ |g k N ω ⬝ᵥ g k N ω - 1|}) (fun N ω hω => ?_)
      (hgdot k 3 (by norm_num))
    have h4 : ¬ (g k N ω ⬝ᵥ g k N ω ≤ 4) := hω
    have h4' : (4 : ℝ) < g k N ω ⬝ᵥ g k N ω := not_le.mp h4
    change (3 : ℝ) ≤ |g k N ω ⬝ᵥ g k N ω - 1|
    rw [abs_of_nonneg (by linarith)]
    linarith
  have hnormvv : ∀ j k : Fin r, Tendsto (fun N => μ N {ω : Ω N |
      WithLp.ofLp (v j N) ⬝ᵥ WithLp.ofLp (v j N) ≤ 4 ∧
        WithLp.ofLp (v k N) ⬝ᵥ WithLp.ofLp (v k N) ≤ 4}ᶜ) atTop (𝓝 0) := by
    intro j k
    have hset : ∀ N : ℕ, {ω : Ω N |
        WithLp.ofLp (v j N) ⬝ᵥ WithLp.ofLp (v j N) ≤ 4 ∧
          WithLp.ofLp (v k N) ⬝ᵥ WithLp.ofLp (v k N) ≤ 4}ᶜ = (∅ : Set (Ω N)) := by
      intro N
      ext ω
      simp [hvdot j N, hvdot k N]
    simp only [hset, measure_empty]
    exact tendsto_const_nhds
  have hnormgg : ∀ j k : Fin r, Tendsto (fun N => μ N {ω : Ω N |
      g j N ω ⬝ᵥ g j N ω ≤ 4 ∧ g k N ω ⬝ᵥ g k N ω ≤ 4}ᶜ) atTop (𝓝 0) := by
    intro j k
    refine tendsto_measure_zero_of_subset (t := fun N =>
      {ω : Ω N | g j N ω ⬝ᵥ g j N ω ≤ 4}ᶜ ∪ {ω : Ω N | g k N ω ⬝ᵥ g k N ω ≤ 4}ᶜ)
      (fun N ω hω => ?_) (tendsto_measure_zero_union (hgbad j) (hgbad k))
    have hnot : ¬ (g j N ω ⬝ᵥ g j N ω ≤ 4 ∧ g k N ω ⬝ᵥ g k N ω ≤ 4) := hω
    by_cases h1 : g j N ω ⬝ᵥ g j N ω ≤ 4
    · exact Set.mem_union_right _ (fun h2 => hnot ⟨h1, h2⟩)
    · exact Set.mem_union_left _ h1
  have hnormvg : ∀ j k : Fin r, Tendsto (fun N => μ N {ω : Ω N |
      WithLp.ofLp (v j N) ⬝ᵥ WithLp.ofLp (v j N) ≤ 4 ∧
        g k N ω ⬝ᵥ g k N ω ≤ 4}ᶜ) atTop (𝓝 0) := by
    intro j k
    refine tendsto_measure_zero_of_subset (fun N ω hω => ?_) (hgbad k)
    have hnot : ¬ (WithLp.ofLp (v j N) ⬝ᵥ WithLp.ofLp (v j N) ≤ 4 ∧
        g k N ω ⬝ᵥ g k N ω ≤ 4) := hω
    refine fun h2 => hnot ⟨?_, h2⟩
    rw [hvdot j N]
    norm_num
  have hzeroL : ∀ x : ℝ, Tendsto
      (fun η : ℝ => (fun _ : ℂ => (0 : ℂ)) ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] (0 : ℝ))
      (𝓝 (((0 : ℝ) : ℂ))) := by
    intro x
    simp only [Complex.ofReal_zero]
    exact tendsto_const_nhds
  -- the six complex limits, per pair
  have hvvC : ∀ (j k : Fin r) (z : ℂ), 0 < z.im → TendstoInProb μ (fun N ω =>
      ‖R4C.cformC (W₀ N ω) z (WithLp.ofLp (v j N)) (WithLp.ofLp (v k N))
        - (if j = k then MP.mC c z else 0)‖) 0 := by
    intro j k z hz
    rcases eq_or_ne j k with rfl | hjk
    · simp only [if_true]
      have h := R2.tendstoInProb_qformC_v (fun N ω => ((ZZ N ω).1 j, (ZZ N ω).2)) hz
        (FormsR.hasLaw_pairRow ZZ hZZ j) (v j) (hnorm j) hd (hR1a z hz)
      simpa only [hW0eq, FormsR.qformC_eq_cformC] using h
    · simp only [if_neg hjk, sub_zero]
      have h := FormsR.tendstoInProb_cformC_cross (μ := μ)
        (Wf := fun N ω => R2.W0 (ZZ N ω).2) (fun N ω => R2.isHermitian_W0 _) hz
        (af := fun N _ => WithLp.ofLp (v j N)) (bf := fun N _ => WithLp.ofLp (v k N))
        (uf := fun N _ => WithLp.ofLp ((Real.sqrt 2)⁻¹ • (v j N + v k N)))
        (L := MP.mC c z) (fun N ω => hvsum j k N)
        (R2.tendstoInProb_qformC_v (fun N ω => ((ZZ N ω).1 j, (ZZ N ω).2)) hz
          (FormsR.hasLaw_pairRow ZZ hZZ j) (v j) (hnorm j) hd (hR1a z hz))
        (R2.tendstoInProb_qformC_v (fun N ω => ((ZZ N ω).1 j, (ZZ N ω).2)) hz
          (FormsR.hasLaw_pairRow ZZ hZZ j) (v k) (hnorm k) hd (hR1a z hz))
        (R2.tendstoInProb_qformC_v (fun N ω => ((ZZ N ω).1 j, (ZZ N ω).2)) hz
          (FormsR.hasLaw_pairRow ZZ hZZ j) (fun N => (Real.sqrt 2)⁻¹ • (v j N + v k N))
          (hup j k hjk) hd (hR1a z hz))
      simpa only [hW0eq] using h
  have hvv2C : ∀ (j k : Fin r) (z : ℂ), 0 < z.im → TendstoInProb μ (fun N ω =>
      ‖R4C.cform2C (W₀ N ω) z (WithLp.ofLp (v j N)) (WithLp.ofLp (v k N))
        - (if j = k then MP.mCDeriv c z else 0)‖) 0 := by
    intro j k z hz
    rcases eq_or_ne j k with rfl | hjk
    · simp only [if_true]
      have h := R2.tendstoInProb_qform2C_v (fun N ω => ((ZZ N ω).1 j, (ZZ N ω).2)) hz
        (FormsR.hasLaw_pairRow ZZ hZZ j) (v j) (hnorm j) hd (hR1b z hz)
      simpa only [hW0eq, FormsR.qform2C_eq_cform2C] using h
    · simp only [if_neg hjk, sub_zero]
      have h := FormsR.tendstoInProb_cform2C_cross (μ := μ)
        (Wf := fun N ω => R2.W0 (ZZ N ω).2) (fun N ω => R2.isHermitian_W0 _) hz
        (af := fun N _ => WithLp.ofLp (v j N)) (bf := fun N _ => WithLp.ofLp (v k N))
        (uf := fun N _ => WithLp.ofLp ((Real.sqrt 2)⁻¹ • (v j N + v k N)))
        (L := MP.mCDeriv c z) (fun N ω => hvsum j k N)
        (R2.tendstoInProb_qform2C_v (fun N ω => ((ZZ N ω).1 j, (ZZ N ω).2)) hz
          (FormsR.hasLaw_pairRow ZZ hZZ j) (v j) (hnorm j) hd (hR1b z hz))
        (R2.tendstoInProb_qform2C_v (fun N ω => ((ZZ N ω).1 j, (ZZ N ω).2)) hz
          (FormsR.hasLaw_pairRow ZZ hZZ j) (v k) (hnorm k) hd (hR1b z hz))
        (R2.tendstoInProb_qform2C_v (fun N ω => ((ZZ N ω).1 j, (ZZ N ω).2)) hz
          (FormsR.hasLaw_pairRow ZZ hZZ j) (fun N => (Real.sqrt 2)⁻¹ • (v j N + v k N))
          (hup j k hjk) hd (hR1b z hz))
      simpa only [hW0eq] using h
  have hggC : ∀ (j k : Fin r) (z : ℂ), 0 < z.im → TendstoInProb μ (fun N ω =>
      ‖R4C.cformC (W₀ N ω) z (g j N ω) (g k N ω) - (if j = k then MP.mC c z else 0)‖) 0 := by
    intro j k z hz
    rcases eq_or_ne j k with rfl | hjk
    · simp only [if_true]
      have h := FormsR.tendstoInProb_qformC_g ZZ hz hZZ hd j (hR1a z hz)
      simpa only [hW0eq, hgeq, FormsR.qformC_eq_cformC] using h
    · simp only [if_neg hjk, sub_zero]
      have h := FormsR.tendstoInProb_cformC_cross (μ := μ)
        (Wf := fun N ω => R2.W0 (ZZ N ω).2) (fun N ω => R2.isHermitian_W0 _) hz
        (af := fun N ω => R2.gOf ((ZZ N ω).1 j)) (bf := fun N ω => R2.gOf ((ZZ N ω).1 k))
        (uf := fun N ω => R2.gOf (fun l => (Real.sqrt 2)⁻¹ *
          ((ZZ N ω).1 j l + 1 * (ZZ N ω).1 k l)))
        (L := MP.mC c z)
        (fun N ω => FormsR.gOf_add_eq_sqrt_two_smul ((ZZ N ω).1 j) ((ZZ N ω).1 k))
        (FormsR.tendstoInProb_qformC_g ZZ hz hZZ hd j (hR1a z hz))
        (FormsR.tendstoInProb_qformC_g ZZ hz hZZ hd k (hR1a z hz))
        (FormsR.tendstoInProb_qformC_gmix ZZ hz hZZ hd hjk (hR1a z hz))
      simpa only [hW0eq, hgeq] using h
  have hgg2C : ∀ (j k : Fin r) (z : ℂ), 0 < z.im → TendstoInProb μ (fun N ω =>
      ‖R4C.cform2C (W₀ N ω) z (g j N ω) (g k N ω)
        - (if j = k then MP.mCDeriv c z else 0)‖) 0 := by
    intro j k z hz
    rcases eq_or_ne j k with rfl | hjk
    · simp only [if_true]
      have h := FormsR.tendstoInProb_qform2C_g ZZ hz hZZ hd j (hR1b z hz)
      simpa only [hW0eq, hgeq, FormsR.qform2C_eq_cform2C] using h
    · simp only [if_neg hjk, sub_zero]
      have h := FormsR.tendstoInProb_cform2C_cross (μ := μ)
        (Wf := fun N ω => R2.W0 (ZZ N ω).2) (fun N ω => R2.isHermitian_W0 _) hz
        (af := fun N ω => R2.gOf ((ZZ N ω).1 j)) (bf := fun N ω => R2.gOf ((ZZ N ω).1 k))
        (uf := fun N ω => R2.gOf (fun l => (Real.sqrt 2)⁻¹ *
          ((ZZ N ω).1 j l + 1 * (ZZ N ω).1 k l)))
        (L := MP.mCDeriv c z)
        (fun N ω => FormsR.gOf_add_eq_sqrt_two_smul ((ZZ N ω).1 j) ((ZZ N ω).1 k))
        (FormsR.tendstoInProb_qform2C_g ZZ hz hZZ hd j (hR1b z hz))
        (FormsR.tendstoInProb_qform2C_g ZZ hz hZZ hd k (hR1b z hz))
        (FormsR.tendstoInProb_qform2C_gmix ZZ hz hZZ hd hjk (hR1b z hz))
      simpa only [hW0eq, hgeq] using h
  have hvgC : ∀ (j k : Fin r) (z : ℂ), 0 < z.im → TendstoInProb μ (fun N ω =>
      ‖R4C.cformC (W₀ N ω) z (WithLp.ofLp (v j N)) (g k N ω)‖) 0 := by
    intro j k z hz
    have h := R2.tendstoInProb_cformC_gOf (fun N ω => ((ZZ N ω).1 k, (ZZ N ω).2)) hz
      (FormsR.hasLaw_pairRow ZZ hZZ k) (v j) (hnorm j) hd
    simpa only [hW0eq, hgeq] using h
  have hvg2C : ∀ (j k : Fin r) (z : ℂ), 0 < z.im → TendstoInProb μ (fun N ω =>
      ‖R4C.cform2C (W₀ N ω) z (WithLp.ofLp (v j N)) (g k N ω)‖) 0 := by
    intro j k z hz
    have h := R2.tendstoInProb_cform2C_gOf (fun N ω => ((ZZ N ω).1 k, (ZZ N ω).2)) hz
      (FormsR.hasLaw_pairRow ZZ hZZ k) (v j) (hnorm j) hd
    simpa only [hW0eq, hgeq] using h
  -- item T, once per field and pair
  refine ⟨hedge, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro j k z hz
    rcases eq_or_ne j k with rfl | hjk
    · simp only [if_true]
      exact T.tendstoInProb_cform_of_complex' hc W₀ hsymm
        (fun N _ => WithLp.ofLp (v j N)) (fun N _ => WithLp.ofLp (v j N)) hedgeC
        (hnormvv j j) (MP.mC c) (MP.m c z)
        (fun z' hz' => by simpa only [if_true] using hvvC j j z' hz')
        hz (MP.tendsto_mC hc hz)
    · simp only [if_neg hjk]
      exact T.tendstoInProb_cform_of_complex' hc W₀ hsymm
        (fun N _ => WithLp.ofLp (v j N)) (fun N _ => WithLp.ofLp (v k N)) hedgeC
        (hnormvv j k) (fun _ => (0 : ℂ)) 0
        (fun z' hz' => by simpa only [if_neg hjk, sub_zero] using hvvC j k z' hz')
        hz (hzeroL z)
  · intro j k z hz
    rcases eq_or_ne j k with rfl | hjk
    · simp only [if_true]
      exact T.tendstoInProb_cform_of_complex' hc W₀ hsymm (g j) (g j) hedgeC
        (hnormgg j j) (MP.mC c) (MP.m c z)
        (fun z' hz' => by simpa only [if_true] using hggC j j z' hz')
        hz (MP.tendsto_mC hc hz)
    · simp only [if_neg hjk]
      exact T.tendstoInProb_cform_of_complex' hc W₀ hsymm (g j) (g k) hedgeC
        (hnormgg j k) (fun _ => (0 : ℂ)) 0
        (fun z' hz' => by simpa only [if_neg hjk, sub_zero] using hggC j k z' hz')
        hz (hzeroL z)
  · intro j k z hz
    exact T.tendstoInProb_cform_of_complex' hc W₀ hsymm
      (fun N _ => WithLp.ofLp (v j N)) (g k) hedgeC (hnormvg j k) (fun _ => (0 : ℂ)) 0
      (fun z' hz' => by simpa only [sub_zero] using hvgC j k z' hz') hz (hzeroL z)
  · intro j k z hz
    rcases eq_or_ne j k with rfl | hjk
    · simp only [if_true]
      exact T.tendstoInProb_cform2_of_complex' hc W₀ hsymm
        (fun N _ => WithLp.ofLp (v j N)) (fun N _ => WithLp.ofLp (v j N)) hedgeC
        (hnormvv j j) (MP.mCDeriv c) (MP.mDeriv c z)
        (fun z' hz' => by simpa only [if_true] using hvv2C j j z' hz')
        hz (MP.tendsto_mCDeriv hc hz)
    · simp only [if_neg hjk]
      exact T.tendstoInProb_cform2_of_complex' hc W₀ hsymm
        (fun N _ => WithLp.ofLp (v j N)) (fun N _ => WithLp.ofLp (v k N)) hedgeC
        (hnormvv j k) (fun _ => (0 : ℂ)) 0
        (fun z' hz' => by simpa only [if_neg hjk, sub_zero] using hvv2C j k z' hz')
        hz (hzeroL z)
  · intro j k z hz
    rcases eq_or_ne j k with rfl | hjk
    · simp only [if_true]
      exact T.tendstoInProb_cform2_of_complex' hc W₀ hsymm (g j) (g j) hedgeC
        (hnormgg j j) (MP.mCDeriv c) (MP.mDeriv c z)
        (fun z' hz' => by simpa only [if_true] using hgg2C j j z' hz')
        hz (MP.tendsto_mCDeriv hc hz)
    · simp only [if_neg hjk]
      exact T.tendstoInProb_cform2_of_complex' hc W₀ hsymm (g j) (g k) hedgeC
        (hnormgg j k) (fun _ => (0 : ℂ)) 0
        (fun z' hz' => by simpa only [if_neg hjk, sub_zero] using hgg2C j k z' hz')
        hz (hzeroL z)
  · intro j k z hz
    exact T.tendstoInProb_cform2_of_complex' hc W₀ hsymm
      (fun N _ => WithLp.ofLp (v j N)) (g k) hedgeC (hnormvg j k) (fun _ => (0 : ℂ)) 0
      (fun z' hz' => by simpa only [sub_zero] using hvg2C j k z' hz') hz (hzeroL z)

end Transfer

/-! ### 6. Real bilinearity helpers and the derived `Q`-column limits -/

namespace FormsR

section RealForms

open R4

variable {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ}

theorem cform_add_left (z : ℝ) (x x' y : Fin D → ℝ) :
    cform W z (x + x') y = cform W z x y + cform W z x' y := by
  simp only [cform, add_dotProduct]

theorem cform2_add_left (z : ℝ) (x x' y : Fin D → ℝ) :
    cform2 W z (x + x') y = cform2 W z x y + cform2 W z x' y := by
  simp only [cform2, add_dotProduct]

theorem cform2_add_right' (z : ℝ) (x y y' : Fin D → ℝ) :
    cform2 W z x (y + y') = cform2 W z x y + cform2 W z x y' := by
  simp only [cform2, Matrix.mulVec_add, dotProduct_add]

theorem cform2_smul_right' (z θ : ℝ) (x y : Fin D → ℝ) :
    cform2 W z x (θ • y) = θ * cform2 W z x y := by
  simp only [cform2, Matrix.mulVec_smul, dotProduct_smul, smul_eq_mul]

/-- The real resolvent form is symmetric. -/
theorem cform_comm (hW : W.IsHermitian) (z : ℝ) (x y : Fin D → ℝ) :
    cform W z x y = cform W z y x :=
  dotProduct_mulVec_comm (transpose_resolv hW) x y

/-- The real squared-resolvent form is symmetric. -/
theorem cform2_comm (hW : W.IsHermitian) (z : ℝ) (x y : Fin D → ℝ) :
    cform2 W z x y = cform2 W z y x := by
  refine dotProduct_mulVec_comm ?_ x y
  rw [Matrix.transpose_mul, transpose_resolv hW]

end RealForms

/-- Move a `TendstoInProb` across an equality of limits. -/
theorem tendstoInProb_congr_limit {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {f : ∀ N, Ω N → ℝ} {a b : ℝ} (hab : a = b)
    (h : TendstoInProb μ f a) : TendstoInProb μ f b := hab ▸ h

end FormsR

section Derived

open R4

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {dN : ℕ → ℕ} {r : ℕ} {c : ℝ}
  {W₀ : (N : ℕ) → Ω N → Matrix (Fin (dN N)) (Fin (dN N)) ℝ}
  {hsymm : ∀ N ω, (W₀ N ω).IsHermitian}
  {v : Fin r → (N : ℕ) → EuclideanSpace ℝ (Fin (dN N))}
  {g : Fin r → (N : ℕ) → Ω N → Fin (dN N) → ℝ}

/-- **The plan's matrix limit, entrywise, first order** (plan 1.2): entry `(j, k)` of
`Qᵀ G₀(z) Q` tends to `δ_jk (λ_k + 1) m(z)` at every real `z > bulkEdge c`. Column `j` of
`Q` is `√λ_j v_j + g_j`. -/
theorem ResolventLimitsR.cform_qcol (h : ResolventLimitsR μ W₀ hsymm v g c)
    {lam : Fin r → ℝ} (hlam : ∀ j, 0 ≤ lam j) (j k : Fin r) {z : ℝ}
    (hz : bulkEdge c < z) :
    TendstoInProb μ (fun N ω => cform (W₀ N ω) z
        (Real.sqrt (lam j) • WithLp.ofLp (v j N) + g j N ω)
        (Real.sqrt (lam k) • WithLp.ofLp (v k N) + g k N ω))
      (if j = k then (lam k + 1) * MP.m c z else 0) := by
  have hfun : (fun N (ω : Ω N) => cform (W₀ N ω) z
        (Real.sqrt (lam j) • WithLp.ofLp (v j N) + g j N ω)
        (Real.sqrt (lam k) • WithLp.ofLp (v k N) + g k N ω))
      = fun N ω =>
        Real.sqrt (lam j) * Real.sqrt (lam k) *
            cform (W₀ N ω) z (WithLp.ofLp (v j N)) (WithLp.ofLp (v k N))
          + Real.sqrt (lam j) * cform (W₀ N ω) z (WithLp.ofLp (v j N)) (g k N ω)
          + Real.sqrt (lam k) * cform (W₀ N ω) z (WithLp.ofLp (v k N)) (g j N ω)
          + cform (W₀ N ω) z (g j N ω) (g k N ω) := by
    funext N ω
    simp only [FormsR.cform_add_left, cform_add_right, cform_smul_left, cform_smul_right]
    rw [FormsR.cform_comm (hsymm N ω) z (g j N ω) (WithLp.ofLp (v k N))]
    ring
  have hcomb := ((((h.vv j k z hz).const_mul
      (Real.sqrt (lam j) * Real.sqrt (lam k))).add
      ((h.vg j k z hz).const_mul (Real.sqrt (lam j)))).add
      ((h.vg k j z hz).const_mul (Real.sqrt (lam k)))).add (h.gg j k z hz)
  rw [hfun]
  rcases eq_or_ne j k with rfl | hjk
  · refine FormsR.tendstoInProb_congr_limit ?_ hcomb
    simp only [if_true]
    rw [Real.mul_self_sqrt (hlam j)]
    ring
  · refine FormsR.tendstoInProb_congr_limit ?_ hcomb
    simp only [if_neg hjk]
    ring

/-- **The `G₀²` twin, entrywise** (plan 1.2): entry `(j, k)` of `Qᵀ G₀(z)² Q` tends to
`δ_jk (λ_k + 1) m'(z)`. -/
theorem ResolventLimitsR.cform2_qcol (h : ResolventLimitsR μ W₀ hsymm v g c)
    {lam : Fin r → ℝ} (hlam : ∀ j, 0 ≤ lam j) (j k : Fin r) {z : ℝ}
    (hz : bulkEdge c < z) :
    TendstoInProb μ (fun N ω => cform2 (W₀ N ω) z
        (Real.sqrt (lam j) • WithLp.ofLp (v j N) + g j N ω)
        (Real.sqrt (lam k) • WithLp.ofLp (v k N) + g k N ω))
      (if j = k then (lam k + 1) * MP.mDeriv c z else 0) := by
  have hfun : (fun N (ω : Ω N) => cform2 (W₀ N ω) z
        (Real.sqrt (lam j) • WithLp.ofLp (v j N) + g j N ω)
        (Real.sqrt (lam k) • WithLp.ofLp (v k N) + g k N ω))
      = fun N ω =>
        Real.sqrt (lam j) * Real.sqrt (lam k) *
            cform2 (W₀ N ω) z (WithLp.ofLp (v j N)) (WithLp.ofLp (v k N))
          + Real.sqrt (lam j) * cform2 (W₀ N ω) z (WithLp.ofLp (v j N)) (g k N ω)
          + Real.sqrt (lam k) * cform2 (W₀ N ω) z (WithLp.ofLp (v k N)) (g j N ω)
          + cform2 (W₀ N ω) z (g j N ω) (g k N ω) := by
    funext N ω
    simp only [FormsR.cform2_add_left, FormsR.cform2_add_right', cform2_smul_left,
      FormsR.cform2_smul_right']
    rw [FormsR.cform2_comm (hsymm N ω) z (g j N ω) (WithLp.ofLp (v k N))]
    ring
  have hcomb := ((((h.vv2 j k z hz).const_mul
      (Real.sqrt (lam j) * Real.sqrt (lam k))).add
      ((h.vg2 j k z hz).const_mul (Real.sqrt (lam j)))).add
      ((h.vg2 k j z hz).const_mul (Real.sqrt (lam k)))).add (h.gg2 j k z hz)
  rw [hfun]
  rcases eq_or_ne j k with rfl | hjk
  · refine FormsR.tendstoInProb_congr_limit ?_ hcomb
    simp only [if_true]
    rw [Real.mul_self_sqrt (hlam j)]
    ring
  · refine FormsR.tendstoInProb_congr_limit ?_ hcomb
    simp only [if_neg hjk]
    ring

end Derived

/-! ### 6b. The shared stack interface: the `r`-indexed assembly

The same assembly as section 7, on `RankRStack` (`RankR/RMT/Stack.lean`) instead of on
`UnalignedModel`. Nothing here reads `r_i = 1`: the inputs are the split of section 4b of
`RankR/RMT/Split.lean`, the aspect ratio of the stack, and the Gaussian law of the noise.
Section 7 is this theorem at `UnalignedModel.toStack`. -/

section StackInterface

open R4

namespace RankRStack

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-- Column `k` of `G = Eᵀ U` is `gOf` of row `k` of `Uᵀ Zu` (the bridge from the split
objects to the pair law, `E_transpose_mul`). -/
theorem gOf_row_eq_col (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (ns N)) (Fin r) ℝ) (k : Fin r) :
    R2.gOf ((Uᵀ * s.Zu N ω) k) = fun l => ((s.E N ω)ᵀ * U) l k := by
  funext l
  rw [s.E_transpose_mul N ω U]
  simp only [Matrix.smul_apply, Matrix.transpose_apply, smul_eq_mul, R2.gOf]

/-- **Task U2, the assembly** (plan section 2) on the shared interface: the `r`-indexed
`ResolventLimitsR` for a `RankRStack` with Gaussian noise, on the split of section 4b.

The side condition that remains is `hn : ∀ N, ns N = r + pp N` together with
`hp : ∀ N, 0 < pp N` (that is, `r < ns N` at every `N`); decision D11's tail-shift removal
is task U7's job. -/
theorem resolventLimitsR_of_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (s : RankRStack μ ns d r) (hG : s.GaussianNoise) {c : ℝ} (hc : 0 < c)
    (hdtop : Tendsto d atTop atTop)
    (hns : Tendsto (fun N => (ns N : ℝ) / d N) atTop (𝓝 c))
    {pp : ℕ → ℕ} (hn : ∀ N, ns N = r + pp N) (hp : ∀ N, 0 < pp N) :
    ∃ U : (N : ℕ) → Matrix (Fin (ns N)) (Fin r) ℝ,
      (∀ N, (U N)ᵀ * U N = 1) ∧
      (∀ N, s.signalPart N
        = U N * (s.spikeMat N * Matrix.diagonal fun j => Real.sqrt (s.coreEig j))ᵀ) ∧
      (∀ N ω, s.gram N ω
        = s.rankRW0 N ω (U N) + s.qmatR N ω (U N) * (s.qmatR N ω (U N))ᵀ) ∧
      ResolventLimitsR μ (fun N ω => s.rankRW0 N ω (U N))
        (fun N ω => s.isHermitian_rankRW0 N ω (U N)) s.spikeVec
        (fun k N ω => fun l => ((s.E N ω)ᵀ * U N) l k) c := by
  classical
  choose U B hU hsig hgram hlaw using fun N => s.exists_rankR_split hG N (hn N)
  have hd0 : ∀ N, 0 < d N := s.hd
  have hr : ∀ N, r ≤ ns N := fun N => by
    rw [hn N]
    exact Nat.le_add_right r (pp N)
  have hppN : ∀ N, ns N - r = pp N := fun N => by
    rw [hn N]
    omega
  have hcN : Tendsto (fun N => (pp N : ℝ) / d N) atTop (𝓝 c) := by
    have h := tendsto_sub_ratio (nf := ns) r hr hdtop hns
    simpa only [hppN] using h
  set ZZ : ∀ N, Ω N → Matrix (Fin r) (Fin (d N)) ℝ × Matrix (Fin (pp N)) (Fin (d N)) ℝ :=
    fun N ω => ((U N)ᵀ * s.Zu N ω, B N ω) with hZZdef
  have hZZlaw : ∀ N, HasLaw (ZZ N)
      ((gaussianMatrix r (d N)).prod (gaussianMatrix (pp N) (d N))) (μ N) := hlaw
  have hBsnd : ∀ N, HasLaw (B N) (gaussianMatrix (pp N) (d N)) (μ N) :=
    fun N => measurePreserving_snd.fun_comp_hasLaw (hlaw N)
  have hW0eq : ∀ N ω, R2.W0 (ZZ N ω).2 = s.rankRW0 N ω (U N) := fun N ω =>
    (s.rankRW0_eq_smul_block N ω (hU N) (hsig N) (hgram N ω)).symm
  have hgeq : ∀ (k : Fin r) (N : ℕ) (ω : Ω N),
      R2.gOf ((ZZ N ω).1 k) = fun l => ((s.E N ω)ᵀ * U N) l k :=
    fun k N ω => s.gOf_row_eq_col N ω (U N) k
  have hedge1 : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
        ≤ bulkEdge c + ε}) atTop (𝓝 1) := by
    intro ε hε
    exact R3.tendsto_measure_lamMax_le μ hc B (fun N ω => s.rankRW0 N ω (U N))
      (fun N ω => s.rankRW0_eq_smul_block N ω (hU N) (hsig N) (hgram N ω))
      (fun N ω => s.isHermitian_rankRW0 N ω (U N)) hBsnd hp hd0 hdtop hcN hε
  have hedgeC : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
        ≤ bulkEdge c + ε}ᶜ) atTop (𝓝 0) := fun ε hε =>
    tendsto_measure_compl_zero
      (fun N => (s.measurableSet_lamMax_rankRW0_le N (U N) _).nullMeasurableSet)
      (hedge1 ε hε)
  refine ⟨U, hU, hsig, fun N ω => s.gram_eq_rankRW0_add N ω (hU N) (hsig N), ?_⟩
  exact resolventLimitsR_of_pairLaw hc hdtop hp hcN (fun N ω => s.rankRW0 N ω (U N))
    (fun N ω => s.isHermitian_rankRW0 N ω (U N)) s.spikeVec
    (fun k N ω => fun l => ((s.E N ω)ᵀ * U N) l k) ZZ hZZlaw hW0eq hgeq
    (fun a b N => s.inner_spikeVec N a b) hedge1 hedgeC

/-- The plan's `Qᵀ G₀(z) Q → (Λ + 1) m(z)`, entrywise: `Q = qmatR`, `Λ = diag(coreEig)`. -/
theorem tendstoInProb_cform_qmatR (s : RankRStack μ ns d r)
    {U : (N : ℕ) → Matrix (Fin (ns N)) (Fin r) ℝ} {c : ℝ}
    (h : ResolventLimitsR μ (fun N ω => s.rankRW0 N ω (U N))
      (fun N ω => s.isHermitian_rankRW0 N ω (U N)) s.spikeVec
      (fun k N ω => fun l => ((s.E N ω)ᵀ * U N) l k) c)
    (j k : Fin r) {z : ℝ} (hz : bulkEdge c < z) :
    TendstoInProb μ (fun N ω => R4.cform (s.rankRW0 N ω (U N)) z
        (fun l => s.qmatR N ω (U N) l j) (fun l => s.qmatR N ω (U N) l k))
      (if j = k then (s.coreEig k + 1) * MP.m c z else 0) := by
  have hq := h.cform_qcol (fun a => s.coreEig_nonneg a) j k hz
  have hcol : ∀ (a : Fin r) (N : ℕ) (ω : Ω N), (fun l => s.qmatR N ω (U N) l a)
      = Real.sqrt (s.coreEig a) • WithLp.ofLp (s.spikeVec a N)
        + fun l => ((s.E N ω)ᵀ * U N) l a := by
    intro a N ω
    funext l
    rw [Pi.add_apply, Pi.smul_apply, smul_eq_mul, s.qmatR_apply N ω (U N) l a]
  simp only [hcol]
  exact hq

/-- The `G₀²` twin of `tendstoInProb_cform_qmatR`. -/
theorem tendstoInProb_cform2_qmatR (s : RankRStack μ ns d r)
    {U : (N : ℕ) → Matrix (Fin (ns N)) (Fin r) ℝ} {c : ℝ}
    (h : ResolventLimitsR μ (fun N ω => s.rankRW0 N ω (U N))
      (fun N ω => s.isHermitian_rankRW0 N ω (U N)) s.spikeVec
      (fun k N ω => fun l => ((s.E N ω)ᵀ * U N) l k) c)
    (j k : Fin r) {z : ℝ} (hz : bulkEdge c < z) :
    TendstoInProb μ (fun N ω => R4.cform2 (s.rankRW0 N ω (U N)) z
        (fun l => s.qmatR N ω (U N) l j) (fun l => s.qmatR N ω (U N) l k))
      (if j = k then (s.coreEig k + 1) * MP.mDeriv c z else 0) := by
  have hq := h.cform2_qcol (fun a => s.coreEig_nonneg a) j k hz
  have hcol : ∀ (a : Fin r) (N : ℕ) (ω : Ω N), (fun l => s.qmatR N ω (U N) l a)
      = Real.sqrt (s.coreEig a) • WithLp.ofLp (s.spikeVec a N)
        + fun l => ((s.E N ω)ᵀ * U N) l a := by
    intro a N ω
    funext l
    rw [Pi.add_apply, Pi.smul_apply, smul_eq_mul, s.qmatR_apply N ω (U N) l a]
  simp only [hcol]
  exact hq

end RankRStack

end StackInterface

/-! ### 7. The model: `W₀` of the split, measurability, and the assembly -/

namespace UnalignedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-- The frame complement of the stacked noise, `E⊥ = E_stack - U Uᵀ E_stack`. -/
noncomputable def stackEperp (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ) :
    Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
  m.stackE N ω - U * (Uᵀ * m.stackE N ω)

/-- The rank-`r` Wishart block `W₀ = E⊥ᵀ E⊥` of the split (plan 1.1 item 2). -/
noncomputable def rankRW0 (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ) : Matrix (Fin (d N)) (Fin (d N)) ℝ :=
  (m.stackEperp N ω U)ᵀ * m.stackEperp N ω U

theorem isHermitian_rankRW0 (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ) : (m.rankRW0 N ω U).IsHermitian :=
  isHermitian_transpose_mul_self _

/-! #### The bridge to `RankRStack` (continued from `RankR/RMT/Split.lean`) -/

theorem toStack_stackEperp [NeZero M] (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ) :
    m.toStack.stackEperp N ω U = m.stackEperp N ω U := rfl

theorem toStack_rankRW0 [NeZero M] (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ) :
    m.toStack.rankRW0 N ω U = m.rankRW0 N ω U := rfl

/-- The stack as frame times coefficients plus complement (plan 1.1 item 2), for any `U`
with the factor property. -/
theorem stackX_eq_frame_add (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    {U : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ}
    (hsig : m.signalPart N
      = U * (m.spikeMat N * Matrix.diagonal fun j => Real.sqrt (m.coreEig j))ᵀ) :
    m.stackX N ω = U * (m.qmatR N ω U)ᵀ + m.stackEperp N ω U := by
  have hT : ((m.stackE N ω)ᵀ * U)ᵀ = Uᵀ * m.stackE N ω := by
    rw [Matrix.transpose_mul, Matrix.transpose_transpose]
  rw [m.stackX_eq N ω, hsig, qmatR, Matrix.transpose_add, Matrix.mul_add, hT, stackEperp]
  abel

/-- The Gram matrix of the stack splits as `W₀ + Q Qᵀ` (plan 1.1 item 2). -/
theorem stackGram_eq_rankRW0_add (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    {U : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ} (hU : Uᵀ * U = 1)
    (hsig : m.signalPart N
      = U * (m.spikeMat N * Matrix.diagonal fun j => Real.sqrt (m.coreEig j))ᵀ) :
    m.stackGram N ω = m.rankRW0 N ω U + m.qmatR N ω U * (m.qmatR N ω U)ᵀ := by
  have hFperp : (m.stackEperp N ω U)ᵀ * U = 0 := transpose_perp_mul U (m.stackE N ω) hU
  rw [stackGram, m.stackX_eq_frame_add N ω hsig,
    gram_mul_transpose_add U (m.qmatR N ω U) (m.stackEperp N ω U) hU hFperp]
  rfl

/-- On the split of `exists_rankR_split` the block identity `W₀ = d⁻¹ Bᵀ B` holds
pointwise. -/
theorem rankRW0_eq_smul_block (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    {U : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ} (hU : Uᵀ * U = 1)
    (hsig : m.signalPart N
      = U * (m.spikeMat N * Matrix.diagonal fun j => Real.sqrt (m.coreEig j))ᵀ)
    {pp : ℕ} {Bm : Matrix (Fin pp) (Fin (d N)) ℝ}
    (hgram : m.stackGram N ω
      = ((d N : ℝ))⁻¹ • (Bmᵀ * Bm) + m.qmatR N ω U * (m.qmatR N ω U)ᵀ) :
    m.rankRW0 N ω U = ((d N : ℝ))⁻¹ • (Bmᵀ * Bm) := by
  have h1 := m.stackGram_eq_rankRW0_add N ω hU hsig
  rw [h1] at hgram
  exact add_right_cancel hgram

/-- Column `k` of `G = E_stackᵀ U` is `gOf` of row `k` of `Uᵀ Z_stack` (the bridge from the
split objects to the pair law, `stackE_transpose_mul`). -/
theorem gOf_row_eq_col (m : UnalignedModel μ M n d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ) (k : Fin r) :
    R2.gOf ((Uᵀ * m.stackZu N ω) k) = fun l => ((m.stackE N ω)ᵀ * U) l k := by
  funext l
  rw [m.stackE_transpose_mul N ω U]
  simp only [Matrix.smul_apply, Matrix.transpose_apply, smul_eq_mul, R2.gOf]

/-- The spike directions are orthonormal, in the inner-product form the interface takes. -/
theorem inner_spikeVec (m : UnalignedModel μ M n d r) (N : ℕ) (a b : Fin r) :
    ⟪m.spikeVec a N, m.spikeVec b N⟫_ℝ = if a = b then 1 else 0 := by
  have h := congrFun (congrFun (m.spikeMat_transpose_mul_self N) a) b
  rw [Matrix.mul_apply] at h
  simp only [Matrix.transpose_apply, Matrix.one_apply] at h
  rw [inner_euclidean_eq_dotProduct]
  calc WithLp.ofLp (m.spikeVec a N) ⬝ᵥ WithLp.ofLp (m.spikeVec b N)
      = ∑ l, m.spikeMat N l a * m.spikeMat N l b := by
        simp only [dotProduct, m.spikeMat_apply N]
    _ = if a = b then 1 else 0 := h

/-- `0 ≤ λ_j(C)`: the core matrix is a Gram matrix. -/
theorem coreEig_nonneg (m : UnalignedModel μ M n d r) (j : Fin r) : 0 ≤ m.coreEig j := by
  have hpsd : m.core.PosSemidef := by
    change ((BR (fun i => (m.tbl i).θ) m.R)ᵀ * BR (fun i => (m.tbl i).θ) m.R).PosSemidef
    simpa using Matrix.posSemidef_conjTranspose_mul_self (BR (fun i => (m.tbl i).θ) m.R)
  exact hpsd.eigenvalues_nonneg j

theorem measurable_stackEperp (m : UnalignedModel μ M n d r) (N : ℕ)
    (U : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ) :
    Measurable (fun ω => m.stackEperp N ω U) := by
  have hE : ∀ (q : Fin (∑ i, n i N)) (l : Fin (d N)),
      Measurable fun ω => m.stackE N ω q l := by
    intro q l
    have hZm : Measurable fun ω => (m.tbl (finSigmaFinEquiv.symm q).1).Z N ω
        (finSigmaFinEquiv.symm q).2 l := by
      have h1 := (m.tbl (finSigmaFinEquiv.symm q).1).hZ N
      exact (measurable_pi_apply l).comp ((measurable_pi_apply _).comp h1)
    have h : (fun ω => m.stackE N ω q l)
        = fun ω => (Real.sqrt (d N))⁻¹ *
          (m.tbl (finSigmaFinEquiv.symm q).1).Z N ω (finSigmaFinEquiv.symm q).2 l := rfl
    rw [h]
    exact hZm.const_mul _
  refine measurable_pi_lambda _ fun q => measurable_pi_lambda _ fun l => ?_
  have h : (fun ω => m.stackEperp N ω U q l)
      = fun ω => m.stackE N ω q l
        - ∑ a : Fin r, U q a * ∑ b, U b a * m.stackE N ω b l := by
    funext ω
    simp only [stackEperp, Matrix.sub_apply, Matrix.mul_apply, Matrix.transpose_apply]
  rw [h]
  refine (hE q l).sub ?_
  refine Finset.measurable_sum _ fun a _ => ?_
  exact (Finset.measurable_sum _ fun b _ => (hE b l).const_mul (U b a)).const_mul (U q a)

theorem measurableSet_lamMax_rankRW0_le (m : UnalignedModel μ M n d r) (N : ℕ)
    (U : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ) (x : ℝ) :
    MeasurableSet {ω | lamMax (m.rankRW0 N ω U) (m.isHermitian_rankRW0 N ω U) ≤ x} := by
  have h : {ω | lamMax (m.rankRW0 N ω U) (m.isHermitian_rankRW0 N ω U) ≤ x}
      = (fun ω => gramLamMax (m.stackEperp N ω U)) ⁻¹' Set.Iic x := rfl
  rw [h]
  exact (measurable_gramLamMax.comp (m.measurable_stackEperp N U)) measurableSet_Iic

/-- **Task U2, the assembly** (plan section 2): the `r`-indexed `ResolventLimitsR` for the
stack of an `UnalignedModel` with joint Gaussian noise, on the split of U1.

The side condition that remains is `hn : ∀ N, ∑ i, n i N = r + pp N` together with
`hp : ∀ N, 0 < pp N` (that is, `r < ∑ i, n i N` at every `N`); its removal by the
`TailShift` device (decision D11) is task U7's job. -/
theorem resolventLimitsR_of_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ M n d r) (hG : m.JointGaussianNoise)
    {cc : Fin M → ℝ} (hreg : ∀ i, (m.tbl i).Regime (cc i)) (hc : 0 < ∑ i, cc i)
    {pp : ℕ → ℕ} (hn : ∀ N, ∑ i, n i N = r + pp N) (hp : ∀ N, 0 < pp N) :
    ∃ U : (N : ℕ) → Matrix (Fin (∑ i, n i N)) (Fin r) ℝ,
      (∀ N, (U N)ᵀ * U N = 1) ∧
      (∀ N, m.signalPart N
        = U N * (m.spikeMat N * Matrix.diagonal fun j => Real.sqrt (m.coreEig j))ᵀ) ∧
      (∀ N ω, m.stackGram N ω
        = m.rankRW0 N ω (U N) + m.qmatR N ω (U N) * (m.qmatR N ω (U N))ᵀ) ∧
      ResolventLimitsR μ (fun N ω => m.rankRW0 N ω (U N))
        (fun N ω => m.isHermitian_rankRW0 N ω (U N)) m.spikeVec
        (fun k N ω => fun l => ((m.stackE N ω)ᵀ * U N) l k) (∑ i, cc i) :=
  m.toStack.resolventLimitsR_of_gaussian (m.gaussianNoise_toStack hG) hc
    (hreg ⟨0, NeZero.pos M⟩).2.1 (m.toMultiTable.stack_regime cc hreg).2.2 hn hp

/-- The plan's `Qᵀ G₀(z) Q → (Λ + 1) m(z)`, entrywise on the model: `Q = qmatR`,
`Λ = diag(coreEig)`. -/
theorem tendstoInProb_cform_qmatR (m : UnalignedModel μ M n d r)
    {U : (N : ℕ) → Matrix (Fin (∑ i, n i N)) (Fin r) ℝ} {c : ℝ}
    (h : ResolventLimitsR μ (fun N ω => m.rankRW0 N ω (U N))
      (fun N ω => m.isHermitian_rankRW0 N ω (U N)) m.spikeVec
      (fun k N ω => fun l => ((m.stackE N ω)ᵀ * U N) l k) c)
    (j k : Fin r) {z : ℝ} (hz : bulkEdge c < z) :
    TendstoInProb μ (fun N ω => R4.cform (m.rankRW0 N ω (U N)) z
        (fun l => m.qmatR N ω (U N) l j) (fun l => m.qmatR N ω (U N) l k))
      (if j = k then (m.coreEig k + 1) * MP.m c z else 0) := by
  have hq := h.cform_qcol (fun a => m.coreEig_nonneg a) j k hz
  have hcol : ∀ (a : Fin r) (N : ℕ) (ω : Ω N), (fun l => m.qmatR N ω (U N) l a)
      = Real.sqrt (m.coreEig a) • WithLp.ofLp (m.spikeVec a N)
        + fun l => ((m.stackE N ω)ᵀ * U N) l a := by
    intro a N ω
    funext l
    rw [Pi.add_apply, Pi.smul_apply, smul_eq_mul, m.qmatR_apply N ω (U N) l a]
  simp only [hcol]
  exact hq

/-- The `G₀²` twin of `tendstoInProb_cform_qmatR`. -/
theorem tendstoInProb_cform2_qmatR (m : UnalignedModel μ M n d r)
    {U : (N : ℕ) → Matrix (Fin (∑ i, n i N)) (Fin r) ℝ} {c : ℝ}
    (h : ResolventLimitsR μ (fun N ω => m.rankRW0 N ω (U N))
      (fun N ω => m.isHermitian_rankRW0 N ω (U N)) m.spikeVec
      (fun k N ω => fun l => ((m.stackE N ω)ᵀ * U N) l k) c)
    (j k : Fin r) {z : ℝ} (hz : bulkEdge c < z) :
    TendstoInProb μ (fun N ω => R4.cform2 (m.rankRW0 N ω (U N)) z
        (fun l => m.qmatR N ω (U N) l j) (fun l => m.qmatR N ω (U N) l k))
      (if j = k then (m.coreEig k + 1) * MP.mDeriv c z else 0) := by
  have hq := h.cform2_qcol (fun a => m.coreEig_nonneg a) j k hz
  have hcol : ∀ (a : Fin r) (N : ℕ) (ω : Ω N), (fun l => m.qmatR N ω (U N) l a)
      = Real.sqrt (m.coreEig a) • WithLp.ofLp (m.spikeVec a N)
        + fun l => ((m.stackE N ω)ᵀ * U N) l a := by
    intro a N ω
    funext l
    rw [Pi.add_apply, Pi.smul_apply, smul_eq_mul, m.qmatR_apply N ω (U N) l a]
  simp only [hcol]
  exact hq

end UnalignedModel

end StackedSVD
