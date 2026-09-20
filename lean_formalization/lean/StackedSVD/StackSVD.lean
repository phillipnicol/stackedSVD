/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT
import StackedSVD.RMT.Full
import StackedSVD.Prob.TendstoInProb

/-!
# `prop:stacksvd_general`: unweighted stacked SVD

STATUS 2026-08-31: proved on 2026-08-29 and 2026-08-30. The file has no `sorry` and no row in
`docs/SORRIES.md`; the Gaussian corollary is discharged by `RMT/Full.lean`. The review note
`notes/archive/prop_stacksvd_general.md` is `user OK (bulk)` (decision D19 item 3).

The paper proves `prop:stacksvd_general` from BEJ1129 Theorem 2.3 at `w = 1_M`. This file
takes the other route of the note: with unit weights and independent Gaussian tables the
stacked matrix is again a rank-one spiked table, with signal `‖θ‖₂` and aspect ratio
`‖c‖₁`. So the black box is `SingleTableLaw` of the stack and the content is the stacking
lemma.

## Content

1. `stackSVDLimit`, the closed-form limit, and `betaSq_sqrt`, the scalar identity that links
   it to `betaSq` at `θ = ‖θ‖₂`.
2. `MultiTableModel.stack`, the unit-weight stack as a `SpikedModel`, with rows
   `Fin (∑ i, n i N)` obtained from the block index `(i : Fin M) × Fin (n i N)` by
   `finSigmaFinEquiv`.
3. The stacking lemmas `stack_X`, `stack_regime`, `stack_law`.
4. `prop_stacksvd_general`, the Layer 1 implication.
5. `prop_stacksvd_general_gaussian`, the Layer 2 corollary: the same conclusion from the
   regime and the joint Gaussian law alone.
6. `prop_stacksvd_general_inner` and `prop_stacksvd_general_inner_gaussian`, the paper's own
   form `|⟨v̂_stacksvd, v⟩|²` for any selected unit top eigenvector of the stack Gram matrix.
   The transfer from the projector form reads `SingleTableLaw.topSimple` of the stack
   (`overlap_eq_inner_sq`), so it holds almost everywhere at each `N`, not pointwise.

`prop:stacksvd_general` and `thm_simple_thm1_stacksvd_gaussian` moved to `StackSVD/Main.lean`
(F28, 2026-09-08); the `private` helper `stackSVDLimit_const_eq` that F28 moved with them was
dropped there by F32 (2026-09-08) for the public `Scalars.simple_thm1_stacksvd`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

/-! ### The scalar limit -/

/-- Closed-form limit of unweighted stacksvd (`prop:stacksvd_general`), written in the
paper's data `θ` and `c`. Equality `‖θ‖₂⁴ = ‖c‖₁` gives `0`, as in `betaSq`. -/
noncomputable def stackSVDLimit {M : ℕ} (θ c : Fin M → ℝ) : ℝ :=
  if (∑ i, θ i ^ 2) ^ 2 > ∑ i, c i then
    ((∑ i, θ i ^ 2) ^ 2 - ∑ i, c i) / ((∑ i, θ i ^ 2) * ((∑ i, θ i ^ 2) + 1))
  else 0

/-- `betaSq` at a square root. With `T = ‖θ‖₂²` and `C = ‖c‖₁` the right side is
`stackSVDLimit`. -/
theorem betaSq_sqrt {T C : ℝ} (hT : 0 ≤ T) :
    betaSq (Real.sqrt T) C = if T ^ 2 > C then (T ^ 2 - C) / (T * (T + 1)) else 0 := by
  have h2 : Real.sqrt T ^ 2 = T := Real.sq_sqrt hT
  have h4 : Real.sqrt T ^ 4 = T ^ 2 := by
    rw [show (4 : ℕ) = 2 * 2 from rfl, pow_mul, h2]
  unfold betaSq
  rw [h4, h2]
  split_ifs with h
  · rw [show T ^ 2 + T = T * (T + 1) by ring]
  · rfl


namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

section Stack

-- One table at least. `Fin M` must be nonempty: the stack needs a row and a shared `v`.
variable [NeZero M]

/-! ### The unit-weight stack -/

/-- `‖θ‖₂²`, the total squared signal strength of the `M` tables. -/
noncomputable def stackThetaSq (m : MultiTableModel μ M n d) : ℝ := ∑ i, (m.tbl i).θ ^ 2

/-- `‖θ‖₂`, the signal strength of the unit-weight stack. -/
noncomputable def stackTheta (m : MultiTableModel μ M n d) : ℝ := Real.sqrt m.stackThetaSq

/-- The stacked left singular vector in block form: block `i` is `θ_i u_i / ‖θ‖₂`. -/
noncomputable def stackUSigma (m : MultiTableModel μ M n d) (N : ℕ) :
    ((i : Fin M) × Fin (n i N)) → ℝ :=
  fun p => (m.tbl p.1).θ * (m.tbl p.1).u N p.2 / m.stackTheta

/-- Left singular vector of the stack. When `‖θ‖₂ = 0` the stacked signal vanishes and any
unit vector serves; the first basis vector of the block index is used (choice 3 of the
note). -/
noncomputable def stackU (m : MultiTableModel μ M n d) (N : ℕ) :
    EuclideanSpace ℝ (Fin (∑ i, n i N)) :=
  if m.stackThetaSq = 0 then
    EuclideanSpace.single (finSigmaFinEquiv ⟨(0 : Fin M), ⟨0, (m.tbl 0).hn N⟩⟩) 1
  else WithLp.toLp 2 fun k => m.stackUSigma N (finSigmaFinEquiv.symm k)

/-- Noise of the stack: the `M` noise matrices placed in one column of blocks, rows
reindexed from `(i : Fin M) × Fin (n i N)` to `Fin (∑ i, n i N)`. -/
noncomputable def stackZ (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
  Matrix.reindex finSigmaFinEquiv (Equiv.refl (Fin (d N)))
    (Matrix.of fun (p : (i : Fin M) × Fin (n i N)) k => (m.tbl p.1).Z N ω p.2 k)

/-- The stack has at least one row. -/
theorem stack_row_pos (m : MultiTableModel μ M n d) (N : ℕ) : 0 < ∑ i, n i N := by
  have : Nonempty (Fin M) := ⟨0⟩
  exact Finset.sum_pos (fun i _ => (m.tbl i).hn N) Finset.univ_nonempty

/-- The stacked left singular vector is a unit vector. -/
theorem norm_stackU (m : MultiTableModel μ M n d) (N : ℕ) : ‖m.stackU N‖ = 1 := by
  rw [stackU]
  split_ifs with h
  · simp
  · have hpos : 0 < m.stackThetaSq :=
      lt_of_le_of_ne (Finset.sum_nonneg fun i _ => sq_nonneg _) (Ne.symm h)
    have hθ2 : m.stackTheta ^ 2 = m.stackThetaSq := Real.sq_sqrt hpos.le
    have hu : ∀ i : Fin M, ∑ j, ((m.tbl i).u N j) ^ 2 = 1 := by
      intro i
      have hnorm := EuclideanSpace.norm_sq_eq ((m.tbl i).u N)
      rw [(m.tbl i).hu N] at hnorm
      simpa [Real.norm_eq_abs, sq_abs] using hnorm.symm
    have key : ∀ i : Fin M, ∑ j, (m.stackUSigma N ⟨i, j⟩) ^ 2
        = (m.tbl i).θ ^ 2 / m.stackThetaSq := by
      intro i
      have hterm : ∀ j : Fin (n i N), (m.stackUSigma N ⟨i, j⟩) ^ 2
          = ((m.tbl i).θ ^ 2 / m.stackThetaSq) * ((m.tbl i).u N j) ^ 2 := by
        intro j
        change ((m.tbl i).θ * (m.tbl i).u N j / m.stackTheta) ^ 2 = _
        rw [div_pow, mul_pow, hθ2]
        ring
      rw [Finset.sum_congr rfl fun j _ => hterm j, ← Finset.mul_sum, hu i, mul_one]
    rw [EuclideanSpace.norm_eq, Real.sqrt_eq_one]
    have hpt : ∀ k : Fin (∑ i, n i N),
        ‖(WithLp.toLp 2 fun k => m.stackUSigma N (finSigmaFinEquiv.symm k) :
          EuclideanSpace ℝ (Fin (∑ i, n i N))) k‖ ^ 2
          = (m.stackUSigma N (finSigmaFinEquiv.symm k)) ^ 2 := by
      intro k
      simp [Real.norm_eq_abs, sq_abs]
    rw [Finset.sum_congr rfl fun k _ => hpt k,
      Equiv.sum_comp finSigmaFinEquiv.symm (fun p => (m.stackUSigma N p) ^ 2),
      Fintype.sum_sigma (fun p => (m.stackUSigma N p) ^ 2),
      Finset.sum_congr rfl fun i _ => key i, ← Finset.sum_div]
    exact div_self (ne_of_gt hpos)

set_option linter.unusedSectionVars false in
/-- The stacked noise is measurable. -/
theorem measurable_stackZ (m : MultiTableModel μ M n d) (N : ℕ) :
    Measurable (m.stackZ N) := by
  refine measurable_pi_lambda _ fun r => measurable_pi_lambda _ fun k => ?_
  exact ((measurable_pi_apply k).comp
    ((measurable_pi_apply (finSigmaFinEquiv.symm r).2).comp
      ((m.tbl (finSigmaFinEquiv.symm r).1).hZ N)))

/-- The unit-weight stack `[X_1; ...; X_M]` as one rank-one spiked table, with signal
`‖θ‖₂` and shared right singular vector `v`. -/
noncomputable def stack (m : MultiTableModel μ M n d) :
    SpikedModel μ (fun N => ∑ i, n i N) d where
  θ := m.stackTheta
  u := m.stackU
  v := (m.tbl 0).v
  Z := m.stackZ
  hθ := Real.sqrt_nonneg _
  hn := m.stack_row_pos
  hd := (m.tbl 0).hd
  hu := m.norm_stackU
  hv := (m.tbl 0).hv
  hZ := m.measurable_stackZ

@[simp]
theorem stack_theta (m : MultiTableModel μ M n d) :
    m.stack.θ = Real.sqrt (∑ i, (m.tbl i).θ ^ 2) := rfl

@[simp]
theorem stack_v (m : MultiTableModel μ M n d) (N : ℕ) : m.stack.v N = (m.tbl 0).v N := rfl

private theorem stack_theta_eq (m : MultiTableModel μ M n d) : m.stack.θ = m.stackTheta := rfl

private theorem stack_u_eq (m : MultiTableModel μ M n d) (N : ℕ) :
    m.stack.u N = m.stackU N := rfl

private theorem stack_Z_eq (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) :
    m.stack.Z N ω = m.stackZ N ω := rfl

/-! ### The stacking lemmas -/

/-- Stacking lemma 1: block `i` of the stacked data matrix is `X_i`. The scaling
`θ_i u_i / ‖θ‖₂` of `stackU` is chosen for this. -/
theorem stack_X (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N)
    (i : Fin M) (j : Fin (n i N)) (k : Fin (d N)) :
    m.stack.X N ω (finSigmaFinEquiv ⟨i, j⟩) k = (m.tbl i).X N ω j k := by
  have hv : (m.tbl 0).v N = (m.tbl i).v N := m.hv 0 i N
  have hZ : m.stackZ N ω (finSigmaFinEquiv ⟨i, j⟩) k = (m.tbl i).Z N ω j k := by
    simp [stackZ, Matrix.reindex_apply, Matrix.submatrix_apply]
  have hu : m.stackTheta * m.stackU N (finSigmaFinEquiv ⟨i, j⟩)
      = (m.tbl i).θ * (m.tbl i).u N j := by
    rw [stackU]
    split_ifs with hzero
    · have hθi : (m.tbl i).θ = 0 := by
        have hsq := (Finset.sum_eq_zero_iff_of_nonneg
          (fun l (_ : l ∈ Finset.univ) => sq_nonneg ((m.tbl l).θ))).1 hzero i (Finset.mem_univ i)
        exact (pow_eq_zero_iff (n := 2) (by norm_num)).1 hsq
      have hst : m.stackTheta = 0 := by
        rw [stackTheta, hzero, Real.sqrt_zero]
      rw [hst, hθi]
      ring
    · have hpos : 0 < m.stackThetaSq :=
        lt_of_le_of_ne (Finset.sum_nonneg fun l _ => sq_nonneg _) (Ne.symm hzero)
      have hne : m.stackTheta ≠ 0 := ne_of_gt (Real.sqrt_pos.mpr hpos)
      rw [PiLp.toLp_apply, Equiv.symm_apply_apply, stackUSigma]
      field_simp
  simp only [SpikedModel.X, SpikedModel.E, Matrix.add_apply, Matrix.smul_apply, smul_eq_mul,
    Matrix.vecMulVec_apply, stack_theta_eq, stack_u_eq, stack_Z_eq, stack_v]
  rw [hZ, hv, ← mul_assoc, hu, mul_assoc]

/-- Stacking lemma 2: the proportional regime carries over, with `c_stack = ∑ i, c i`. -/
theorem stack_regime (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (h : ∀ i, (m.tbl i).Regime (c i)) : m.stack.Regime (∑ i, c i) := by
  refine ⟨?_, (h 0).2.1, ?_⟩
  · refine tendsto_atTop_mono (fun N => ?_) (h 0).1
    exact Finset.single_le_sum (f := fun l => n l N) (fun l _ => Nat.zero_le _)
      (Finset.mem_univ 0)
  · have hcast : ∀ N : ℕ, ((∑ i, n i N : ℕ) : ℝ) / (d N : ℝ) = ∑ i, ((n i N : ℝ) / (d N : ℝ)) := by
      intro N
      rw [Nat.cast_sum, Finset.sum_div]
    simp only [hcast]
    exact tendsto_finsetSum _ fun i _ => (h i).2.2

set_option linter.unusedSectionVars false in
/-- The stacking map on plain product types: the family `Zs` becomes one array whose row `r`
is row `(finSigmaFinEquiv.symm r).2` of table `(finSigmaFinEquiv.symm r).1`. Stated on
`(i : Fin M) → Fin (n i N) → Fin (d N) → ℝ` rather than on `Matrix`, because instance search
does not unfold `Matrix` and so does not find the product sigma-algebra there. -/
private theorem measurePreserving_stackPi (N : ℕ) :
    MeasurePreserving
      (fun (Zs : (i : Fin M) → Fin (n i N) → Fin (d N) → ℝ) (r : Fin (∑ i, n i N)) =>
        Zs (finSigmaFinEquiv.symm r).1 (finSigmaFinEquiv.symm r).2)
      (Measure.pi fun i : Fin M =>
        Measure.pi fun _ : Fin (n i N) => Measure.pi fun _ : Fin (d N) => gaussianReal 0 1)
      (Measure.pi fun _ : Fin (∑ i, n i N) =>
        Measure.pi fun _ : Fin (d N) => gaussianReal 0 1) := by
  have hmeas : Measurable
      (fun (Zs : (i : Fin M) → Fin (n i N) → Fin (d N) → ℝ) (r : Fin (∑ i, n i N)) =>
        Zs (finSigmaFinEquiv.symm r).1 (finSigmaFinEquiv.symm r).2) :=
    measurable_pi_lambda _ fun r =>
      (measurable_pi_apply (finSigmaFinEquiv.symm r).2).comp
        (measurable_pi_apply (finSigmaFinEquiv.symm r).1)
  refine ⟨hmeas, ?_⟩
  refine (Measure.pi_eq fun s hs => ?_).symm
  rw [Measure.map_apply hmeas (MeasurableSet.univ_pi hs)]
  have hpre :
      (fun (Zs : (i : Fin M) → Fin (n i N) → Fin (d N) → ℝ) (r : Fin (∑ i, n i N)) =>
          Zs (finSigmaFinEquiv.symm r).1 (finSigmaFinEquiv.symm r).2) ⁻¹' Set.univ.pi s
        = Set.univ.pi fun i : Fin M =>
            Set.univ.pi fun j : Fin (n i N) => s (finSigmaFinEquiv ⟨i, j⟩) := by
    ext Zs
    simp only [Set.mem_preimage, Set.mem_univ_pi]
    constructor
    · intro hZ i j
      have h1 := hZ (finSigmaFinEquiv ⟨i, j⟩)
      rw [Equiv.symm_apply_apply] at h1
      exact h1
    · intro hZ r
      obtain ⟨p, rfl⟩ : ∃ p, finSigmaFinEquiv p = r :=
        ⟨finSigmaFinEquiv.symm r, Equiv.apply_symm_apply _ _⟩
      obtain ⟨i, j⟩ := p
      rw [Equiv.symm_apply_apply]
      exact hZ i j
  have hRHS : ∏ r : Fin (∑ i, n i N),
        (Measure.pi fun _ : Fin (d N) => gaussianReal 0 1) (s r)
      = ∏ i : Fin M, ∏ j : Fin (n i N),
          (Measure.pi fun _ : Fin (d N) => gaussianReal 0 1) (s (finSigmaFinEquiv ⟨i, j⟩)) := by
    rw [← Equiv.prod_comp finSigmaFinEquiv
      (fun r => (Measure.pi fun _ : Fin (d N) => gaussianReal 0 1) (s r))]
    exact Fintype.prod_sigma _
  rw [hpre, hRHS]
  simp only [Measure.pi_pi]

set_option linter.unusedSectionVars false in
/-- The same statement in the `Matrix` types of `Defs.lean`. -/
private theorem measurePreserving_stackFun (N : ℕ) :
    MeasurePreserving
      (fun (Zs : (i : Fin M) → Matrix (Fin (n i N)) (Fin (d N)) ℝ) =>
        Matrix.reindex finSigmaFinEquiv (Equiv.refl (Fin (d N)))
          (Matrix.of fun (p : (i : Fin M) × Fin (n i N)) k => Zs p.1 p.2 k))
      (Measure.pi fun i : Fin M => gaussianMatrix (n i N) (d N))
      (gaussianMatrix (∑ i, n i N) (d N)) :=
  measurePreserving_stackPi N

/-- Stacking lemma 3: Gaussian noise carries over. The joint law is needed; per-table
`GaussianNoise` does not give independence across tables (audit section 3.8). -/
theorem stack_law (m : MultiTableModel μ M n d) (h : m.JointGaussianNoise) :
    m.stack.GaussianNoise := by
  intro N
  exact (measurePreserving_stackFun N).fun_comp_hasLaw (h N)

end Stack

end MultiTableModel

end StackedSVD
