/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.R3minus
import StackedSVD.RMT.R2

/-!
# Item R6': the subcritical align field

Plan: `notes/archive/rmt_Sym.md`, section "Second use: the subcritical `align` (item R6')". For
`θ⁴ ≤ c` the paper's `β²` is `0`, and this file proves the matching statement

`TendstoInProb μ (fun N ω => overlap (X N ω) (v N)) (betaSq θ c)`.

The route, with R0's split `Xᵀ X = W₀ + q qᵀ`, `q = θ v + g`, `λ := gramLamMax X`:

1. **Decomposition.** `v = a q + r` with `a = ⟪v, q⟫/(q ⬝ q)` and `r ⊥ q`, so
   `overlap X v ≤ 2 a² overlap X q + 2 overlap X r` (`overlap_le_decomp`).
2. **The `q` factor.** On the almost sure simple event there are two cases. If
   `lamMax W₀ < λ` then `λ` is the secular root, `q ⬝ G₀(λ) q = -1`, and
   `R4.topProj_norm_sq` gives `overlap X q = 1 / qform2 W₀ λ q`; `R4.qform2_antitoneOn`
   then bounds it by `1 / qform2 W₀ z₀ q` for every `z₀ ≥ λ`. If `λ = lamMax W₀` the top
   eigenvector of `W₀` is orthogonal to `q` and spans the top space of `Xᵀ X`, so
   `overlap X q = 0` (`overlap_q_le`).
3. **The `r` factor.** Conditionally on `g` the Gram matrix is `q qᵀ + d⁻¹ Bᵀ B` with `B` an
   independent Gaussian block (R0 `exists_block_hasLaw`), which is the Gram matrix of
   `Ymat q d^{-1/2} B`, the matrix with first row `q` and remaining rows `d^{-1/2} B`. Right
   rotations that fix `q` move the test direction and not the law, so item Sym's Bessel
   argument runs with `q` in place of `v` and gives `E[overlap X r | g] ≤ 1/(d-1)`. Fubini on
   R0's product law and Markov turn this into `overlap X r → 0` in probability.
4. **The limit.** `q ⬝ q → θ² + 1` (R2) and `qform2 W₀ z₀ q → (θ²+1) m'(z₀)` ((H2) at the one
   real point `z₀`), so on an event of probability tending to `1`
   `2 a² overlap X q ≤ 8 / ((θ²+1)² m'(z₀))`. `MP.mDeriv_tendsto_atTop` makes this smaller
   than any target once `z₀` is close enough to the bulk edge, so `z₀` is chosen **inside**
   the proof, before the limit in `N` (the bound is vacuous at a fixed `z₀`).

STATUS: see `notes/archive/agent_reports/proof_r6prime.md`.
-/

open Filter Topology MeasureTheory ProbabilityTheory
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD
namespace R6

/-! ### 1. Elementary facts about `overlap` -/

section Elementary

variable {n n' d : ℕ}

/-- `overlap` is homogeneous of degree two in the test direction. -/
theorem overlap_smul (X : Matrix (Fin n) (Fin d) ℝ) (a : ℝ) (w : EuclideanSpace ℝ (Fin d)) :
    overlap X (a • w) = a ^ 2 * overlap X w := by
  change ‖topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) (a • w)‖ ^ 2
      = a ^ 2 * ‖topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) w‖ ^ 2
  rw [map_smul, norm_smul, mul_pow, Real.norm_eq_abs, sq_abs]

theorem overlap_zero (X : Matrix (Fin n) (Fin d) ℝ) :
    overlap X (0 : EuclideanSpace ℝ (Fin d)) = 0 := by
  change ‖topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) 0‖ ^ 2 = 0
  rw [map_zero, norm_zero]
  norm_num

/-- The triangle inequality for the top projector, in squared form. -/
theorem overlap_smul_add_le (X : Matrix (Fin n) (Fin d) ℝ) (a : ℝ)
    (y r : EuclideanSpace ℝ (Fin d)) :
    overlap X (a • y + r) ≤ 2 * (a ^ 2 * overlap X y) + 2 * overlap X r := by
  set P := topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) with hP
  have h1 : ‖P (a • y + r)‖ ≤ |a| * ‖P y‖ + ‖P r‖ := by
    rw [map_add, map_smul]
    refine (norm_add_le _ _).trans ?_
    rw [norm_smul, Real.norm_eq_abs]
  change ‖P (a • y + r)‖ ^ 2 ≤ 2 * (a ^ 2 * ‖P y‖ ^ 2) + 2 * ‖P r‖ ^ 2
  calc ‖P (a • y + r)‖ ^ 2 ≤ (|a| * ‖P y‖ + ‖P r‖) ^ 2 :=
        pow_le_pow_left₀ (norm_nonneg _) h1 2
    _ ≤ 2 * (a ^ 2 * ‖P y‖ ^ 2) + 2 * ‖P r‖ ^ 2 := by
        nlinarith [sq_nonneg (|a| * ‖P y‖ - ‖P r‖), sq_abs a]

end Elementary

/-! ### 2. The orthogonal decomposition of `v` along `q` -/

section Decomposition

variable {d : ℕ}

/-- The coefficient of `q` in `v`. Junk value `0` at `q = 0`. -/
noncomputable def aOf (v q : Fin d → ℝ) : ℝ := (v ⬝ᵥ q) / (q ⬝ᵥ q)

/-- The part of `v` orthogonal to `q`. -/
noncomputable def rvecOf (v q : Fin d → ℝ) : Fin d → ℝ := v - aOf v q • q

/-- Nonnegativity of the dot product with itself. -/
theorem dotProduct_self_nonneg' (q : Fin d → ℝ) : 0 ≤ q ⬝ᵥ q :=
  Finset.sum_nonneg fun _ _ => mul_self_nonneg _

theorem aOf_mul (v q : Fin d → ℝ) : aOf v q * (q ⬝ᵥ q) = v ⬝ᵥ q := by
  rcases eq_or_ne (q ⬝ᵥ q) 0 with h0 | h0
  · have hq : q = 0 := dotProduct_self_eq_zero.mp h0
    subst hq
    simp [aOf]
  · rw [aOf, div_mul_cancel₀ _ h0]

theorem dotProduct_rvecOf (v q : Fin d → ℝ) : rvecOf v q ⬝ᵥ q = 0 := by
  rw [rvecOf, sub_dotProduct, smul_dotProduct, smul_eq_mul, aOf_mul, sub_self]

theorem eq_smul_add_rvecOf (v q : Fin d → ℝ) : v = aOf v q • q + rvecOf v q := by
  rw [rvecOf]
  abel

theorem dotProduct_rvecOf_self_le (v q : Fin d → ℝ) :
    rvecOf v q ⬝ᵥ rvecOf v q ≤ v ⬝ᵥ v := by
  have hexp : rvecOf v q ⬝ᵥ rvecOf v q
      = v ⬝ᵥ v - 2 * (aOf v q * (v ⬝ᵥ q)) + aOf v q ^ 2 * (q ⬝ᵥ q) := by
    rw [rvecOf]
    simp only [sub_dotProduct, dotProduct_sub, smul_dotProduct, dotProduct_smul, smul_eq_mul]
    rw [dotProduct_comm q v]
    ring
  have hac : aOf v q * (v ⬝ᵥ q) = aOf v q ^ 2 * (q ⬝ᵥ q) := by
    rw [← aOf_mul v q]; ring
  rw [hexp, hac]
  nlinarith [mul_nonneg (sq_nonneg (aOf v q)) (dotProduct_self_nonneg' q)]

theorem aOf_sq_mul_le (v q : Fin d → ℝ) : aOf v q ^ 2 * (q ⬝ᵥ q) ≤ v ⬝ᵥ v := by
  have h := dotProduct_rvecOf_self_le v q
  have hnn : 0 ≤ rvecOf v q ⬝ᵥ rvecOf v q := dotProduct_self_nonneg' _
  have hexp : rvecOf v q ⬝ᵥ rvecOf v q
      = v ⬝ᵥ v - aOf v q ^ 2 * (q ⬝ᵥ q) := by
    have hexp' : rvecOf v q ⬝ᵥ rvecOf v q
        = v ⬝ᵥ v - 2 * (aOf v q * (v ⬝ᵥ q)) + aOf v q ^ 2 * (q ⬝ᵥ q) := by
      rw [rvecOf]
      simp only [sub_dotProduct, dotProduct_sub, smul_dotProduct, dotProduct_smul, smul_eq_mul]
      rw [dotProduct_comm q v]
      ring
    have hac : aOf v q * (v ⬝ᵥ q) = aOf v q ^ 2 * (q ⬝ᵥ q) := by
      rw [← aOf_mul v q]; ring
    rw [hexp', hac]; ring
  linarith [hexp ▸ hnn]

theorem toLp_eq_smul_add (v q : Fin d → ℝ) :
    (WithLp.toLp 2 v : EuclideanSpace ℝ (Fin d))
      = aOf v q • WithLp.toLp 2 q + WithLp.toLp 2 (rvecOf v q) := by
  apply WithLp.ofLp_injective
  change v = aOf v q • q + rvecOf v q
  exact eq_smul_add_rvecOf v q

theorem norm_toLp_rvecOf_le (v q : Fin d → ℝ) (hv : v ⬝ᵥ v ≤ 1) :
    ‖(WithLp.toLp 2 (rvecOf v q) : EuclideanSpace ℝ (Fin d))‖ ≤ 1 := by
  have h1 : ‖(WithLp.toLp 2 (rvecOf v q) : EuclideanSpace ℝ (Fin d))‖ ^ 2
      = rvecOf v q ⬝ᵥ rvecOf v q := by
    rw [← real_inner_self_eq_norm_sq]
    exact inner_euclidean_eq_dotProduct _ _
  nlinarith [norm_nonneg (WithLp.toLp 2 (rvecOf v q) : EuclideanSpace ℝ (Fin d)),
    dotProduct_rvecOf_self_le v q, h1]

/-- **Step 1.** The decomposition bound. -/
theorem overlap_le_decomp {n : ℕ} (X : Matrix (Fin n) (Fin d) ℝ) (v q : Fin d → ℝ) :
    overlap X (WithLp.toLp 2 v)
      ≤ 2 * (aOf v q ^ 2 * overlap X (WithLp.toLp 2 q))
        + 2 * overlap X (WithLp.toLp 2 (rvecOf v q)) := by
  rw [toLp_eq_smul_add v q]
  exact overlap_smul_add_le X _ _ _

end Decomposition

/-! ### 3. The deterministic bound on `overlap X q` -/

section QBound

variable {n d : ℕ}

/-- **Step 2, the secular case.** At a secular root above `lamMax W₀` the overlap with `q`
is the inverse of the squared resolvent form: the numerator of `R4.topProj_norm_sq` is
`(q ⬝ G₀ q)² = 1`. -/
theorem overlap_q_eq_inv {W₀ : Matrix (Fin d) (Fin d) ℝ} (hW₀ : W₀.IsHermitian)
    {q : Fin d → ℝ} (hq : q ≠ 0) {lam : ℝ} (hlam : lamMax W₀ hW₀ < lam)
    (e : R4.secular W₀ q lam = 0) (hA : (W₀ + Matrix.vecMulVec q q).IsHermitian) :
    ‖topProj (W₀ + Matrix.vecMulVec q q) hA (WithLp.toLp 2 q)‖ ^ 2
      = 1 / R4.qform2 W₀ lam q := by
  rw [R4.topProj_norm_sq hW₀ hq hlam e hA (WithLp.toLp 2 q)]
  have h2 : q ⬝ᵥ (R4.resolv W₀ lam *ᵥ q) = -1 := by
    have he : 1 + R4.qform W₀ lam q = 0 := e
    have : R4.qform W₀ lam q = q ⬝ᵥ (R4.resolv W₀ lam *ᵥ q) := rfl
    linarith [this ▸ he]
  change (q ⬝ᵥ (R4.resolv W₀ lam *ᵥ q)) ^ 2 / _ = _
  rw [h2]
  norm_num
  rfl

/-- **Step 2.** On the simple event the overlap with `q` is at most the inverse of the
squared resolvent form at any `z₀` above the top eigenvalue of `Xᵀ X`. The case
`gramLamMax X = lamMax W₀` is the one the roadmap left open: there the top eigenvector of
`W₀` is orthogonal to `q` and spans the top space, so the overlap is `0`. -/
theorem overlap_q_le (X : Matrix (Fin n) (Fin d) ℝ) {W₀ : Matrix (Fin d) (Fin d) ℝ}
    (hW₀ : W₀.IsHermitian) {q : Fin d → ℝ}
    (hgram : Xᵀ * X = W₀ + Matrix.vecMulVec q q) (hd : 0 < d)
    (hsimple : TopSimple (Xᵀ * X) (isHermitian_transpose_mul_self X))
    {z₀ : ℝ} (hz₀ : gramLamMax X ≤ z₀) (hpos : 0 < R4.qform2 W₀ z₀ q) :
    overlap X (WithLp.toLp 2 q) ≤ 1 / R4.qform2 W₀ z₀ q := by
  have hA : (W₀ + Matrix.vecMulVec q q).IsHermitian := by
    rw [← hgram]; exact isHermitian_transpose_mul_self X
  have hgl : gramLamMax X = lamMax (W₀ + Matrix.vecMulVec q q) hA :=
    lamMax_congr hgram (isHermitian_transpose_mul_self X) hA
  have hle : lamMax W₀ hW₀ ≤ gramLamMax X := by
    rw [hgl]; exact R4.lamMax_le_lamMax_vecMulVec hW₀ hA
  rcases lt_or_eq_of_le hle with hlt | heq
  · -- the secular case
    have hq : q ≠ 0 := by
      intro h0
      have hz : Matrix.vecMulVec q q = 0 := by
        ext i j; rw [h0]; simp
      have : gramLamMax X = lamMax W₀ hW₀ := by
        rw [hgl]
        exact (lamMax_congr (by rw [hz, add_zero]) hA hW₀)
      linarith
    have hspec : gramLamMax X ∈ spectrum ℝ (toOp (W₀ + Matrix.vecMulVec q q)) := by
      obtain ⟨j, hj⟩ := R4.exists_eigenvalues_eq_lamMax hA hd
      rw [R4.spectrum_toOp, hgl]
      exact hj ▸ hA.eigenvalues_mem_spectrum_real j
    have hsec : R4.secular W₀ q (gramLamMax X) = 0 :=
      (R4.secular_eq_zero_iff hW₀ hlt).mpr hspec
    have hov : overlap X (WithLp.toLp 2 q) = 1 / R4.qform2 W₀ (gramLamMax X) q := by
      change ‖topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) (WithLp.toLp 2 q)‖ ^ 2 = _
      rw [topProj_congr hgram (isHermitian_transpose_mul_self X) hA (WithLp.toLp 2 q)]
      exact overlap_q_eq_inv hW₀ hq hlt hsec hA
    rw [hov]
    have hanti := R4.qform2_antitoneOn hW₀ q (Set.mem_Ioi.2 hlt)
      (Set.mem_Ioi.2 (lt_of_lt_of_le hlt hz₀)) hz₀
    exact one_div_le_one_div_of_le hpos hanti
  · -- the degenerate case: the top eigenvector of `W₀` is orthogonal to `q`
    obtain ⟨x, hx0, hxe⟩ := lamMax_mem_eigSet hW₀ hd
    have hxnorm : ‖x‖ ≠ 0 := norm_ne_zero_iff.mpr hx0
    set φ : EuclideanSpace ℝ (Fin d) := ‖x‖⁻¹ • x with hφdef
    have hφn : ‖φ‖ = 1 := by
      rw [hφdef, norm_smul, norm_inv, norm_norm, inv_mul_cancel₀ hxnorm]
    have hxm : W₀ *ᵥ WithLp.ofLp x = lamMax W₀ hW₀ • WithLp.ofLp x :=
      (R4.mem_eigenspace_iff' W₀ (lamMax W₀ hW₀) x).mp
        (Module.End.mem_eigenspace_iff.mpr hxe)
    have hφe : W₀ *ᵥ WithLp.ofLp φ = lamMax W₀ hW₀ • WithLp.ofLp φ := by
      change W₀ *ᵥ (‖x‖⁻¹ • WithLp.ofLp x) = lamMax W₀ hW₀ • (‖x‖⁻¹ • WithLp.ofLp x)
      rw [Matrix.mulVec_smul, hxm, smul_comm]
    have hφφ : WithLp.ofLp φ ⬝ᵥ WithLp.ofLp φ = 1 := R2.dotProduct_ofLp_self hφn
    have hquad : WithLp.ofLp φ ⬝ᵥ ((W₀ + Matrix.vecMulVec q q) *ᵥ WithLp.ofLp φ)
        = lamMax W₀ hW₀ + (q ⬝ᵥ WithLp.ofLp φ) ^ 2 := by
      rw [Matrix.add_mulVec, dotProduct_add, hφe, R4.vecMulVec_mulVec]
      simp only [dotProduct_smul, smul_eq_mul]
      rw [hφφ, dotProduct_comm (WithLp.ofLp φ) q]
      ring
    have hle2 : WithLp.ofLp φ ⬝ᵥ ((W₀ + Matrix.vecMulVec q q) *ᵥ WithLp.ofLp φ)
        ≤ lamMax (W₀ + Matrix.vecMulVec q q) hA * (WithLp.ofLp φ ⬝ᵥ WithLp.ofLp φ) :=
      R4.dotProduct_mulVec_le_lamMax hA _
    have hqφ : q ⬝ᵥ WithLp.ofLp φ = 0 := by
      rw [hquad, hφφ, mul_one, ← hgl, ← heq] at hle2
      nlinarith [sq_nonneg (q ⬝ᵥ WithLp.ofLp φ)]
    have hAφ : (Xᵀ * X) *ᵥ WithLp.ofLp φ = gramLamMax X • WithLp.ofLp φ := by
      rw [hgram, Matrix.add_mulVec, hφe, R4.vecMulVec_mulVec, hqφ, zero_smul, add_zero, heq]
    have hmem : φ ∈ topSpace (Xᵀ * X) (isHermitian_transpose_mul_self X) := by
      change φ ∈ specSpace (Xᵀ * X) {lamMax (Xᵀ * X) (isHermitian_transpose_mul_self X)}
      rw [specSpace]
      simp only [Set.mem_singleton_iff, iSup_iSup_eq_left]
      exact (R4.mem_eigenspace_iff' _ _ _).mpr hAφ
    have hzero : overlap X (WithLp.toLp 2 q) = 0 := by
      rw [overlap_eq_inner_sq X (WithLp.toLp 2 q) hsimple hmem hφn,
        inner_euclidean_eq_dotProduct]
      change (WithLp.ofLp φ ⬝ᵥ q) ^ 2 = 0
      rw [dotProduct_comm, hqφ]
      norm_num
    rw [hzero]
    exact div_nonneg zero_le_one hpos.le

end QBound

/-! ### 4. The conditional model: rank one plus an independent Gaussian block -/

section Ymat

variable {p d : ℕ}

/-- The `(p+1) × d` matrix with first row `y` and remaining rows `t B i`. Its Gram matrix is
`y yᵀ + t² Bᵀ B`, which is `Xᵀ X` of item R0 at `y = q` and `t = d^{-1/2}`. -/
noncomputable def Ymat (y : Fin d → ℝ) (t : ℝ) (B : Matrix (Fin p) (Fin d) ℝ) :
    Matrix (Fin (p + 1)) (Fin d) ℝ :=
  Matrix.of fun i => Fin.cases y (fun i' => t • B i') i

theorem Ymat_zero (y : Fin d → ℝ) (t : ℝ) (B : Matrix (Fin p) (Fin d) ℝ) :
    Ymat y t B 0 = y := rfl

theorem Ymat_succ (y : Fin d → ℝ) (t : ℝ) (B : Matrix (Fin p) (Fin d) ℝ) (i : Fin p) :
    Ymat y t B i.succ = t • B i := rfl

theorem gram_Ymat (y : Fin d → ℝ) (t : ℝ) (B : Matrix (Fin p) (Fin d) ℝ) :
    (Ymat y t B)ᵀ * Ymat y t B = Matrix.vecMulVec y y + t ^ 2 • (Bᵀ * B) := by
  ext j j'
  change ∑ i : Fin (p + 1), Ymat y t B i j * Ymat y t B i j'
      = y j * y j' + t ^ 2 * ∑ i : Fin p, B i j * B i j'
  have h1 : ∀ i : Fin p, Ymat y t B i.succ j * Ymat y t B i.succ j'
      = t ^ 2 * (B i j * B i j') := by
    intro i
    change (t * B i j) * (t * B i j') = t ^ 2 * (B i j * B i j')
    ring
  rw [Fin.sum_univ_succ, Ymat_zero]
  simp only [h1]
  rw [← Finset.mul_sum]

theorem Ymat_mul (y : Fin d → ℝ) (t : ℝ) (B : Matrix (Fin p) (Fin d) ℝ)
    {O : Matrix (Fin d) (Fin d) ℝ} (hO : Oᵀ *ᵥ y = y) :
    Ymat y t B * O = Ymat y t (B * O) := by
  ext i j
  induction i using Fin.cases with
  | zero =>
      change ∑ k, y k * O k j = y j
      have h := congrFun hO j
      have h' : ∑ k, O k j * y k = y j := by
        simpa [Matrix.mulVec, dotProduct, Matrix.transpose_apply] using h
      rw [← h']
      exact Finset.sum_congr rfl fun k _ => mul_comm _ _
  | succ i' =>
      change ∑ k, (t * B i' k) * O k j = t * ∑ k, B i' k * O k j
      rw [Finset.mul_sum]
      exact Finset.sum_congr rfl fun k _ => by ring

theorem measurable_entry_of (r : Fin p) (j : Fin d) :
    Measurable fun B : Matrix (Fin p) (Fin d) ℝ => B r j := by
  have h1 : Measurable fun B : Matrix (Fin p) (Fin d) ℝ => B r := measurable_pi_apply r
  exact (measurable_pi_apply j).comp h1

theorem measurable_Ymat (y : Fin d → ℝ) (t : ℝ) :
    Measurable (fun B : Matrix (Fin p) (Fin d) ℝ => Ymat y t B) := by
  refine measurable_pi_lambda _ fun i => ?_
  induction i using Fin.cases with
  | zero => exact measurable_const
  | succ i' =>
      refine measurable_pi_lambda _ fun j => ?_
      change Measurable fun B : Matrix (Fin p) (Fin d) ℝ => t * B i' j
      exact (measurable_entry_of i' j).const_mul t

theorem measurable_overlap_Ymat (y : Fin d → ℝ) (t : ℝ) (w : EuclideanSpace ℝ (Fin d)) :
    Measurable (fun B : Matrix (Fin p) (Fin d) ℝ =>
      ENNReal.ofReal (overlap (Ymat y t B) w)) :=
  ENNReal.measurable_ofReal.comp ((measurable_overlap w).comp (measurable_Ymat y t))

/-- **Step 3, exchangeability.** Conditionally on `q` the mean overlap is the same for every
unit direction orthogonal to `q`: the Householder reflection with axis `w - w'` fixes `q` and
`B ↦ B O` preserves the Gaussian law. -/
theorem lintegral_overlap_Ymat_eq (y : Fin d → ℝ) (t : ℝ)
    {w w' : EuclideanSpace ℝ (Fin d)} (hw : ‖w‖ = 1) (hw' : ‖w'‖ = 1)
    (hwy : WithLp.ofLp w ⬝ᵥ y = 0) (hw'y : WithLp.ofLp w' ⬝ᵥ y = 0) :
    ∫⁻ B, ENNReal.ofReal (overlap (Ymat y t B) w) ∂(gaussianMatrix p d)
      = ∫⁻ B, ENNReal.ofReal (overlap (Ymat y t B) w') ∂(gaussianMatrix p d) := by
  by_cases hww : w = w'
  · rw [hww]
  have hdot : ∀ x z : EuclideanSpace ℝ (Fin d),
      WithLp.ofLp x ⬝ᵥ WithLp.ofLp z = ⟪x, z⟫_ℝ := fun x z =>
    (inner_euclidean_eq_dotProduct x z).symm
  have hww1 : WithLp.ofLp w ⬝ᵥ WithLp.ofLp w = 1 := R2.dotProduct_ofLp_self hw
  have hww'1 : WithLp.ofLp w' ⬝ᵥ WithLp.ofLp w' = 1 := R2.dotProduct_ofLp_self hw'
  have hne : WithLp.ofLp w ≠ WithLp.ofLp w' := fun h => hww (by
    have := congrArg (WithLp.toLp 2) h
    simpa using this)
  set a : Fin d → ℝ := WithLp.ofLp w - WithLp.ofLp w' with ha
  set O : Matrix (Fin d) (Fin d) ℝ := householder a with hOdef
  have ha0 : a ⬝ᵥ a ≠ 0 := fun h => (sub_ne_zero.mpr hne) (dotProduct_self_eq_zero.mp h)
  have hO : Oᵀ * O = 1 := householder_orth ha0
  have hay : a ⬝ᵥ y = 0 := by rw [ha, sub_dotProduct, hwy, hw'y, sub_zero]
  have hOy : Oᵀ *ᵥ y = y := by
    rw [hOdef, householder_transpose]
    exact householder_apply_of_orth hay
  have hOw : O *ᵥ WithLp.ofLp w = WithLp.ofLp w' :=
    householder_sub_apply hww1 hww'1 hne
  have key : ∀ B : Matrix (Fin p) (Fin d) ℝ,
      overlap (Ymat y t (B * O)) w = overlap (Ymat y t B) w' := by
    intro B
    rw [← Ymat_mul y t B hOy, overlap_mul_right _ hO, hOw]
  calc ∫⁻ B, ENNReal.ofReal (overlap (Ymat y t B) w) ∂(gaussianMatrix p d)
      = ∫⁻ B, ENNReal.ofReal (overlap (Ymat y t (B * O)) w) ∂(gaussianMatrix p d) :=
        ((measurePreserving_mul_right hO p).lintegral_comp
          (measurable_overlap_Ymat y t w)).symm
    _ = ∫⁻ B, ENNReal.ofReal (overlap (Ymat y t B) w') ∂(gaussianMatrix p d) := by
        simp only [key]

/-- **Step 3, the mean bound.** Bessel over an orthonormal basis of `q^⊥`, which has at least
`d - 1` elements, together with exchangeability. -/
theorem lintegral_overlap_Ymat_le (hd : 2 ≤ d) (y : Fin d → ℝ) (t : ℝ)
    (hsimple : ∀ᵐ B ∂(gaussianMatrix p d),
      TopSimple ((Ymat y t B)ᵀ * Ymat y t B) (isHermitian_transpose_mul_self (Ymat y t B)))
    {w : EuclideanSpace ℝ (Fin d)} (hw : ‖w‖ = 1) (hwy : WithLp.ofLp w ⬝ᵥ y = 0) :
    ∫⁻ B, ENNReal.ofReal (overlap (Ymat y t B) w) ∂(gaussianMatrix p d)
      ≤ ENNReal.ofReal (1 / ((d : ℝ) - 1)) := by
  classical
  set v : EuclideanSpace ℝ (Fin d) := WithLp.toLp 2 y with hvdef
  have hofv : WithLp.ofLp v = y := rfl
  have hrk : d - 1 ≤ Module.finrank ℝ ((ℝ ∙ v)ᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin d))) := by
    have h2 := Submodule.finrank_add_finrank_orthogonal (𝕜 := ℝ) (K := (ℝ ∙ v))
    have h1 : Module.finrank ℝ (ℝ ∙ v) ≤ 1 := by
      rcases eq_or_ne v 0 with h0 | h0
      · rw [h0, Submodule.span_zero_singleton]
        simp
      · rw [finrank_span_singleton h0]
    rw [finrank_euclideanSpace_fin] at h2
    omega
  set b := stdOrthonormalBasis ℝ ((ℝ ∙ v)ᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin d))) with hb
  set e : Fin (Module.finrank ℝ ((ℝ ∙ v)ᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin d)))) →
      EuclideanSpace ℝ (Fin d) := fun k => (b k : EuclideanSpace ℝ (Fin d)) with he
  have hon : Orthonormal ℝ e :=
    (((ℝ ∙ v)ᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin d))).subtypeₗᵢ.orthonormal_comp_iff
      (v := fun k => b k)).mpr b.orthonormal
  have henorm : ∀ k, ‖e k‖ = 1 := fun k => hon.1 k
  have hev : ∀ k, WithLp.ofLp (e k) ⬝ᵥ y = 0 := by
    intro k
    have hmem : (b k : EuclideanSpace ℝ (Fin d)) ∈ (ℝ ∙ v)ᗮ := (b k).2
    rw [Submodule.mem_orthogonal_singleton_iff_inner_right] at hmem
    have hi := (inner_euclidean_eq_dotProduct (e k) v).symm
    rw [hofv] at hi
    rw [hi, real_inner_comm]
    exact hmem
  set I : ℝ≥0∞ := ∫⁻ B, ENNReal.ofReal (overlap (Ymat y t B) w) ∂(gaussianMatrix p d) with hI
  have hsame : ∀ k, ∫⁻ B, ENNReal.ofReal (overlap (Ymat y t B) (e k)) ∂(gaussianMatrix p d)
      = I := fun k =>
    (lintegral_overlap_Ymat_eq y t hw (henorm k) hwy (hev k)).symm
  have hae : ∀ᵐ B ∂(gaussianMatrix p d),
      ENNReal.ofReal (∑ k, overlap (Ymat y t B) (e k)) ≤ 1 := by
    filter_upwards [hsimple] with B hB
    exact ENNReal.ofReal_le_one.mpr (sum_overlap_le_one _ hB hon)
  have hsum : ∑ k, ∫⁻ B, ENNReal.ofReal (overlap (Ymat y t B) (e k)) ∂(gaussianMatrix p d)
      ≤ 1 := by
    rw [← lintegral_finsetSum _ fun k _ => measurable_overlap_Ymat y t (e k)]
    calc ∫⁻ B, ∑ k, ENNReal.ofReal (overlap (Ymat y t B) (e k)) ∂(gaussianMatrix p d)
        = ∫⁻ B, ENNReal.ofReal (∑ k, overlap (Ymat y t B) (e k)) ∂(gaussianMatrix p d) := by
          refine lintegral_congr fun B => ?_
          rw [ENNReal.ofReal_sum_of_nonneg fun k _ => overlap_nonneg _ _]
      _ ≤ ∫⁻ _, (1 : ℝ≥0∞) ∂(gaussianMatrix p d) := lintegral_mono_ae hae
      _ = 1 := by simp
  simp only [hsame, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at hsum
  have hpos : (0 : ℝ) < (d : ℝ) - 1 := by
    have : (2 : ℝ) ≤ (d : ℝ) := by exact_mod_cast hd
    linarith
  have hcard : ENNReal.ofReal ((d : ℝ) - 1)
      ≤ (Module.finrank ℝ ((ℝ ∙ v)ᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin d))) : ℝ≥0∞) := by
    have h1 : ((d : ℝ) - 1)
        ≤ (Module.finrank ℝ ((ℝ ∙ v)ᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin d))) : ℝ) := by
      have hle : ((d - 1 : ℕ) : ℝ)
          ≤ (Module.finrank ℝ ((ℝ ∙ v)ᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin d))) : ℝ) := by
        exact_mod_cast hrk
      have hd1 : (1 : ℕ) ≤ d := by omega
      have : ((d - 1 : ℕ) : ℝ) = (d : ℝ) - 1 := by
        push_cast [Nat.cast_sub hd1]; ring
      linarith [this ▸ hle]
    calc ENNReal.ofReal ((d : ℝ) - 1)
        ≤ ENNReal.ofReal
            ((Module.finrank ℝ ((ℝ ∙ v)ᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin d))) : ℝ)) :=
          ENNReal.ofReal_le_ofReal h1
      _ = _ := by rw [ENNReal.ofReal_natCast]
  have hstep : ENNReal.ofReal ((d : ℝ) - 1) * I ≤ 1 := le_trans (by gcongr) hsum
  have hne0 : ENNReal.ofReal ((d : ℝ) - 1) ≠ 0 := by
    simp only [ne_eq, ENNReal.ofReal_eq_zero, not_le]
    exact hpos
  have hdiv : I ≤ 1 / ENNReal.ofReal ((d : ℝ) - 1) := by
    rw [ENNReal.le_div_iff_mul_le (Or.inl hne0) (Or.inl ENNReal.ofReal_ne_top), mul_comm]
    exact hstep
  refine hdiv.trans (le_of_eq ?_)
  rw [ENNReal.ofReal_div_of_pos hpos, ENNReal.ofReal_one]

/-- The same bound for a direction of norm at most one, by homogeneity. -/
theorem lintegral_overlap_Ymat_le' (hd : 2 ≤ d) (y : Fin d → ℝ) (t : ℝ)
    (hsimple : ∀ᵐ B ∂(gaussianMatrix p d),
      TopSimple ((Ymat y t B)ᵀ * Ymat y t B) (isHermitian_transpose_mul_self (Ymat y t B)))
    {w : EuclideanSpace ℝ (Fin d)} (hw : ‖w‖ ≤ 1) (hwy : WithLp.ofLp w ⬝ᵥ y = 0) :
    ∫⁻ B, ENNReal.ofReal (overlap (Ymat y t B) w) ∂(gaussianMatrix p d)
      ≤ ENNReal.ofReal (1 / ((d : ℝ) - 1)) := by
  rcases eq_or_ne w 0 with h0 | h0
  · subst h0
    simp only [overlap_zero, ENNReal.ofReal_zero, lintegral_zero]
    exact zero_le
  · have hn0 : ‖w‖ ≠ 0 := norm_ne_zero_iff.mpr h0
    set u : EuclideanSpace ℝ (Fin d) := ‖w‖⁻¹ • w with hu
    have hun : ‖u‖ = 1 := by
      rw [hu, norm_smul, norm_inv, norm_norm, inv_mul_cancel₀ hn0]
    have huy : WithLp.ofLp u ⬝ᵥ y = 0 := by
      change (‖w‖⁻¹ • WithLp.ofLp w) ⬝ᵥ y = 0
      rw [smul_dotProduct, hwy, smul_eq_mul, mul_zero]
    have hwu : w = ‖w‖ • u := by
      rw [hu, smul_smul, mul_inv_cancel₀ hn0, one_smul]
    have hmono : ∀ B : Matrix (Fin p) (Fin d) ℝ,
        ENNReal.ofReal (overlap (Ymat y t B) w)
          ≤ ENNReal.ofReal (overlap (Ymat y t B) u) := by
      intro B
      refine ENNReal.ofReal_le_ofReal ?_
      rw [hwu, overlap_smul]
      have h1 : ‖w‖ ^ 2 ≤ 1 := by nlinarith [norm_nonneg w]
      nlinarith [overlap_nonneg (Ymat y t B) u, sq_nonneg ‖w‖]
    calc ∫⁻ B, ENNReal.ofReal (overlap (Ymat y t B) w) ∂(gaussianMatrix p d)
        ≤ ∫⁻ B, ENNReal.ofReal (overlap (Ymat y t B) u) ∂(gaussianMatrix p d) :=
          lintegral_mono hmono
      _ ≤ ENNReal.ofReal (1 / ((d : ℝ) - 1)) :=
          lintegral_overlap_Ymat_le hd y t hsimple hun huy

end Ymat

end R6

/-! ### 5. The model: the `r` direction, its mean overlap, and the align field -/

namespace SpikedModel

open R6

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- The random test direction of the decomposition: the part of `v` orthogonal to `q`. -/
noncomputable def rdir (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) : Fin (d N) → ℝ :=
  R6.rvecOf (WithLp.ofLp (m.v N)) (m.qvec N ω)

/-- `g` in the coordinates R2 uses. -/
theorem gvec_eq_gOf (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) :
    m.gvec N ω = R2.gOf ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)) := by
  rw [m.gvec_eq_smul N ω]
  rfl

theorem qvec_eq (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) :
    m.qvec N ω = m.θ • WithLp.ofLp (m.v N) + R2.gOf ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)) := by
  change m.θ • WithLp.ofLp (m.v N) + m.gvec N ω = _
  rw [gvec_eq_gOf]

/-- The Gram matrix of `X`, written through R0's block. -/
theorem gram_eq_Ymat (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) {p : ℕ}
    {B : Matrix (Fin p) (Fin (d N)) ℝ} (hB : m.W0 N ω = ((d N : ℝ))⁻¹ • (Bᵀ * B)) :
    (m.X N ω)ᵀ * m.X N ω
      = (R6.Ymat (m.qvec N ω) (Real.sqrt (d N))⁻¹ B)ᵀ
          * R6.Ymat (m.qvec N ω) (Real.sqrt (d N))⁻¹ B := by
  have hsq : ((Real.sqrt (d N))⁻¹ : ℝ) ^ 2 = ((d N : ℝ))⁻¹ := by
    rw [sq]
    exact R2.sqrt_inv_mul_sqrt_inv (d N)
  rw [m.gram_eq N ω, R6.gram_Ymat, hsq, hB, add_comm]

/-! #### Measurability of the random objects -/

theorem measurable_gvec (m : SpikedModel μ n d) (N : ℕ) : Measurable (m.gvec N) := by
  have hZ : ∀ (r : Fin (n N)) (j : Fin (d N)), Measurable fun ω => m.Z N ω r j := by
    intro r j
    have h1 : Measurable fun ω => m.Z N ω r := (measurable_pi_apply r).comp (m.hZ N)
    exact (measurable_pi_apply j).comp h1
  have hE : ∀ (r : Fin (n N)) (j : Fin (d N)), Measurable fun ω => m.E N ω r j := by
    intro r j
    have h : (fun ω => m.E N ω r j) = fun ω => (Real.sqrt (d N))⁻¹ * m.Z N ω r j := rfl
    rw [h]
    exact (hZ r j).const_mul _
  refine measurable_pi_lambda _ fun j => ?_
  have h : (fun ω => m.gvec N ω j)
      = fun ω => ∑ k : Fin (n N), m.E N ω k j * WithLp.ofLp (m.u N) k := rfl
  rw [h]
  exact Finset.measurable_sum _ fun k _ => (hE k j).mul_const _

theorem measurable_qvec (m : SpikedModel μ n d) (N : ℕ) : Measurable (m.qvec N) := by
  refine measurable_pi_lambda _ fun j => ?_
  have h : (fun ω => m.qvec N ω j)
      = fun ω => m.θ * WithLp.ofLp (m.v N) j + m.gvec N ω j := rfl
  rw [h]
  exact measurable_const.add ((measurable_pi_apply j).comp (m.measurable_gvec N))

theorem measurable_rvecOf {D : ℕ} (v : Fin D → ℝ) :
    Measurable fun q : Fin D → ℝ => R6.rvecOf v q := by
  have hd1 : Measurable fun q : Fin D → ℝ => v ⬝ᵥ q :=
    Finset.measurable_sum _ fun j _ => (measurable_pi_apply j).const_mul _
  have hd2 : Measurable fun q : Fin D → ℝ => q ⬝ᵥ q :=
    Finset.measurable_sum _ fun j _ => (measurable_pi_apply j).mul (measurable_pi_apply j)
  have hc : Measurable fun q : Fin D → ℝ => R6.aOf v q := hd1.div hd2
  refine measurable_pi_lambda _ fun j => ?_
  change Measurable fun q : Fin D → ℝ => v j - R6.aOf v q * q j
  exact measurable_const.sub (hc.mul (measurable_pi_apply j))

theorem measurable_rdir (m : SpikedModel μ n d) (N : ℕ) :
    Measurable fun ω => (WithLp.toLp 2 (m.rdir N ω) : EuclideanSpace ℝ (Fin (d N))) :=
  (WithLp.measurable_toLp 2 (Fin (d N) → ℝ)).comp
    ((measurable_rvecOf (WithLp.ofLp (m.v N))).comp (m.measurable_qvec N))

theorem measurable_overlap_rdir (m : SpikedModel μ n d) (N : ℕ) :
    Measurable fun ω => overlap (m.X N ω) (WithLp.toLp 2 (m.rdir N ω)) := by
  have h1 : Measurable fun ω : Ω N =>
      (m.X N ω, (WithLp.toLp 2 (m.rdir N ω) : EuclideanSpace ℝ (Fin (d N)))) :=
    (m.measurable_X N).prodMk (m.measurable_rdir N)
  have hf := measurable_overlap₂.comp h1
  exact hf

/-! #### The conditional mean bound and Markov -/

/-- **Step 3, the model.** The mean overlap with the random direction `r` is at most
`1/(d-1)`. The proof conditions on `g` through R0's product law: for a fixed first
coordinate the Gram matrix is `q qᵀ + d⁻¹ Bᵀ B`, item Sym's argument applies with `q` in
place of `v`, and Fubini puts the pieces together. -/
theorem lintegral_overlap_rdir_le (m : SpikedModel μ n d) (hG : m.GaussianNoise) (N : ℕ)
    (hd2 : 2 ≤ d N) :
    ∫⁻ ω, ENNReal.ofReal (overlap (m.X N ω) (WithLp.toLp 2 (m.rdir N ω))) ∂(μ N)
      ≤ ENNReal.ofReal (1 / ((d N : ℝ) - 1)) := by
  classical
  obtain ⟨p, hp⟩ := Nat.exists_eq_succ_of_ne_zero (m.hn N).ne'
  obtain ⟨B, hW, hlaw⟩ := m.exists_block_hasLaw N hG hp
  set t : ℝ := (Real.sqrt (d N))⁻¹ with htdef
  set qOf : (Fin (d N) → ℝ) → (Fin (d N) → ℝ) :=
    fun a => m.θ • WithLp.ofLp (m.v N) + R2.gOf a with hqOfdef
  set F : R2.NoiseSpace p (d N) → ℝ≥0∞ := fun ξ =>
    ENNReal.ofReal (overlap (R6.Ymat (qOf ξ.1) t ξ.2)
      (WithLp.toLp 2 (R6.rvecOf (WithLp.ofLp (m.v N)) (qOf ξ.1)))) with hFdef
  have hgOfm : Measurable fun a : Fin (d N) → ℝ => R2.gOf a := by
    refine measurable_pi_lambda _ fun j => ?_
    change Measurable fun a : Fin (d N) → ℝ => (Real.sqrt (d N))⁻¹ * a j
    exact (measurable_pi_apply j).const_mul _
  have hqm : Measurable qOf := by
    refine measurable_pi_lambda _ fun j => ?_
    change Measurable fun a : Fin (d N) → ℝ =>
      m.θ * WithLp.ofLp (m.v N) j + R2.gOf a j
    exact measurable_const.add ((measurable_pi_apply j).comp hgOfm)
  have hYm : Measurable fun ξ : R2.NoiseSpace p (d N) => R6.Ymat (qOf ξ.1) t ξ.2 := by
    refine measurable_pi_lambda _ fun i => ?_
    induction i using Fin.cases with
    | zero => exact hqm.comp measurable_fst
    | succ i' =>
        refine measurable_pi_lambda _ fun j => ?_
        change Measurable fun ξ : R2.NoiseSpace p (d N) => t * ξ.2 i' j
        exact ((R6.measurable_entry_of i' j).comp measurable_snd).const_mul _
  have hFm : Measurable F := by
    have hrm : Measurable fun ξ : R2.NoiseSpace p (d N) =>
        (WithLp.toLp 2 (R6.rvecOf (WithLp.ofLp (m.v N)) (qOf ξ.1))
          : EuclideanSpace ℝ (Fin (d N))) :=
      (WithLp.measurable_toLp 2 (Fin (d N) → ℝ)).comp
        ((measurable_rvecOf (WithLp.ofLp (m.v N))).comp (hqm.comp measurable_fst))
    have h1 : Measurable fun ξ : R2.NoiseSpace p (d N) =>
        (R6.Ymat (qOf ξ.1) t ξ.2,
          (WithLp.toLp 2 (R6.rvecOf (WithLp.ofLp (m.v N)) (qOf ξ.1))
            : EuclideanSpace ℝ (Fin (d N)))) := hYm.prodMk hrm
    have hf := measurable_overlap₂.comp h1
    have hf2 := ENNReal.measurable_ofReal.comp hf
    exact hf2
  have hq : ∀ ω, qOf ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)) = m.qvec N ω := fun ω =>
    (m.qvec_eq N ω).symm
  have hpt : ∀ ω, ENNReal.ofReal (overlap (m.X N ω) (WithLp.toLp 2 (m.rdir N ω)))
      = F ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N), B ω) := by
    intro ω
    change ENNReal.ofReal (overlap (m.X N ω)
        (WithLp.toLp 2 (R6.rvecOf (WithLp.ofLp (m.v N)) (m.qvec N ω))))
      = ENNReal.ofReal (overlap (R6.Ymat (qOf ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N))) t (B ω))
        (WithLp.toLp 2 (R6.rvecOf (WithLp.ofLp (m.v N))
          (qOf ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N))))))
    rw [hq ω]
    congr 1
    exact overlap_congr_gram _ _ (m.gram_eq_Ymat N ω (hW ω)) _
  -- almost sure simplicity, conditionally on the first coordinate
  have hpred : Measurable fun ξ : R2.NoiseSpace p (d N) =>
      TopSimple ((R6.Ymat (qOf ξ.1) t ξ.2)ᵀ * R6.Ymat (qOf ξ.1) t ξ.2)
        (isHermitian_transpose_mul_self (R6.Ymat (qOf ξ.1) t ξ.2)) := by
    rw [← measurableSet_setOfPred]
    have hset : {ξ : R2.NoiseSpace p (d N) |
        TopSimple ((R6.Ymat (qOf ξ.1) t ξ.2)ᵀ * R6.Ymat (qOf ξ.1) t ξ.2)
          (isHermitian_transpose_mul_self (R6.Ymat (qOf ξ.1) t ξ.2))}
        = (fun ξ : R2.NoiseSpace p (d N) => R6.Ymat (qOf ξ.1) t ξ.2) ⁻¹'
          {Y : Matrix (Fin (p + 1)) (Fin (d N)) ℝ |
            TopSimple (Yᵀ * Y) (isHermitian_transpose_mul_self Y)} := rfl
    rw [hset]
    exact hYm measurableSet_topSimple
  have hprod : ∀ᵐ ξ ∂((Measure.pi fun _ : Fin (d N) => gaussianReal 0 1).prod
      (gaussianMatrix p (d N))),
      TopSimple ((R6.Ymat (qOf ξ.1) t ξ.2)ᵀ * R6.Ymat (qOf ξ.1) t ξ.2)
        (isHermitian_transpose_mul_self (R6.Ymat (qOf ξ.1) t ξ.2)) := by
    rw [← hlaw.ae_iff hpred]
    filter_upwards [singleTableLaw_topSimple_of_gaussian m hG N] with ω hω
    refine topSimple_congr ?_ _ (isHermitian_transpose_mul_self (m.X N ω)) hω
    rw [hq ω]
    exact (m.gram_eq_Ymat N ω (hW ω)).symm
  have hcond := Measure.ae_ae_of_ae_prod hprod
  have hinner : ∀ᵐ a ∂(R2.piGauss (d N)),
      ∫⁻ Bm, F (a, Bm) ∂(gaussianMatrix p (d N))
        ≤ ENNReal.ofReal (1 / ((d N : ℝ) - 1)) := by
    filter_upwards [hcond] with a ha
    have hvv : WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) = 1 :=
      R2.dotProduct_ofLp_self (m.hv N)
    have hFa : ∀ Bm : Matrix (Fin p) (Fin (d N)) ℝ, F (a, Bm)
        = ENNReal.ofReal (overlap (R6.Ymat (qOf a) t Bm)
            (WithLp.toLp 2 (R6.rvecOf (WithLp.ofLp (m.v N)) (qOf a)))) := fun _ => rfl
    simp only [hFa]
    refine R6.lintegral_overlap_Ymat_le' hd2 (qOf a) t ha
      (R6.norm_toLp_rvecOf_le _ _ (le_of_eq hvv)) ?_
    change R6.rvecOf (WithLp.ofLp (m.v N)) (qOf a) ⬝ᵥ qOf a = 0
    exact R6.dotProduct_rvecOf _ _
  calc ∫⁻ ω, ENNReal.ofReal (overlap (m.X N ω) (WithLp.toLp 2 (m.rdir N ω))) ∂(μ N)
      = ∫⁻ ω, F ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N), B ω) ∂(μ N) := lintegral_congr hpt
    _ = ∫⁻ ξ, F ξ ∂((R2.piGauss (d N)).prod (gaussianMatrix p (d N))) :=
        hlaw.lintegral_comp hFm.aemeasurable
    _ = ∫⁻ a, ∫⁻ Bm, F (a, Bm) ∂(gaussianMatrix p (d N)) ∂(R2.piGauss (d N)) :=
        lintegral_prod F hFm.aemeasurable
    _ ≤ ∫⁻ _, ENNReal.ofReal (1 / ((d N : ℝ) - 1)) ∂(R2.piGauss (d N)) :=
        lintegral_mono_ae hinner
    _ = ENNReal.ofReal (1 / ((d N : ℝ) - 1)) := by simp

/-- Markov turns the mean bound into a tail bound. -/
theorem measure_overlap_rdir_ge_le (m : SpikedModel μ n d) (hG : m.GaussianNoise) (N : ℕ)
    (hd2 : 2 ≤ d N) {ε : ℝ} (hε : 0 < ε) :
    μ N {ω | ε ≤ overlap (m.X N ω) (WithLp.toLp 2 (m.rdir N ω))}
      ≤ ENNReal.ofReal (1 / (ε * ((d N : ℝ) - 1))) := by
  have hmeas : AEMeasurable
      (fun ω => ENNReal.ofReal (overlap (m.X N ω) (WithLp.toLp 2 (m.rdir N ω)))) (μ N) :=
    (ENNReal.measurable_ofReal.comp (m.measurable_overlap_rdir N)).aemeasurable
  have hset : {ω | ε ≤ overlap (m.X N ω) (WithLp.toLp 2 (m.rdir N ω))}
      = {ω | ENNReal.ofReal ε
          ≤ ENNReal.ofReal (overlap (m.X N ω) (WithLp.toLp 2 (m.rdir N ω)))} := by
    ext ω
    simp [ENNReal.ofReal_le_ofReal_iff (overlap_nonneg _ _)]
  have hεne : ENNReal.ofReal ε ≠ 0 := by simp [ENNReal.ofReal_eq_zero, not_le, hε]
  have hmk := meas_ge_le_lintegral_div hmeas hεne ENNReal.ofReal_ne_top
  rw [hset]
  refine hmk.trans ?_
  refine (ENNReal.div_le_div_right (m.lintegral_overlap_rdir_le hG N hd2) _).trans
    (le_of_eq ?_)
  rw [← ENNReal.ofReal_div_of_pos hε]
  congr 1
  field_simp

/-- **Step 3, the conclusion.** The overlap with the random direction `r` vanishes in
probability. -/
theorem tendstoInProb_overlap_rdir (m : SpikedModel μ n d) (hG : m.GaussianNoise)
    (hd : Tendsto d atTop atTop) :
    TendstoInProb μ (fun N ω => overlap (m.X N ω) (WithLp.toLp 2 (m.rdir N ω))) 0 := by
  intro ε hε
  have hset : ∀ N : ℕ,
      {ω : Ω N | ε ≤ |overlap (m.X N ω) (WithLp.toLp 2 (m.rdir N ω)) - 0|}
        = {ω : Ω N | ε ≤ overlap (m.X N ω) (WithLp.toLp 2 (m.rdir N ω))} := by
    intro N
    ext ω
    simp [abs_of_nonneg (overlap_nonneg _ _)]
  simp only [hset]
  have hle : ∀ᶠ N in atTop,
      μ N {ω | ε ≤ overlap (m.X N ω) (WithLp.toLp 2 (m.rdir N ω))}
        ≤ ENNReal.ofReal (1 / (ε * ((d N : ℝ) - 1))) := by
    filter_upwards [hd (Filter.eventually_ge_atTop 2)] with N hN
    exact m.measure_overlap_rdir_ge_le hG N hN hε
  have hreal : Tendsto (fun N => 1 / (ε * ((d N : ℝ) - 1))) atTop (𝓝 0) := by
    have h1 : Tendsto (fun N => ε * ((d N : ℝ) - 1)) atTop atTop := by
      refine Filter.Tendsto.const_mul_atTop hε ?_
      exact (tendsto_natCast_atTop_atTop.comp hd).atTop_add tendsto_const_nhds
    exact Filter.Tendsto.congr (fun N => (one_div _).symm) h1.inv_tendsto_atTop
  have htend : Tendsto (fun N => ENNReal.ofReal (1 / (ε * ((d N : ℝ) - 1)))) atTop (𝓝 0) := by
    have h3 := (ENNReal.continuous_ofReal.tendsto 0).comp hreal
    rw [ENNReal.ofReal_zero] at h3
    exact h3
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds htend ?_ hle
  exact Filter.Eventually.of_forall fun _ => by simp

/-! #### The two scalar limits at the model -/

/-- `g ⬝ g → 1`: R2's second moment through R0's block law. -/
theorem tendstoInProb_dotProduct_gvec (m : SpikedModel μ n d) (hG : m.GaussianNoise)
    (hd : Tendsto d atTop atTop) :
    TendstoInProb μ (fun N ω => m.gvec N ω ⬝ᵥ m.gvec N ω) 1 := by
  classical
  have hpsucc : ∀ N, n N = (n N - 1) + 1 := fun N =>
    (Nat.succ_pred_eq_of_pos (m.hn N)).symm
  choose B hW hlaw using fun N => m.exists_block_hasLaw N hG (hpsucc N)
  have h := R2.tendstoInProb_dotProduct_gOf (μ := μ) (pN := fun N => n N - 1) (dN := d)
    (fun N ω => ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N), B N ω)) hlaw hd
  have heq : ∀ (N : ℕ) (ω : Ω N),
      R2.gOf ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)) = m.gvec N ω := fun N ω =>
    (m.gvec_eq_gOf N ω).symm
  simpa only [heq] using h

/-- `v ⬝ g → 0`: the linear form in the Gaussian vector has second moment `1/d`. -/
theorem tendstoInProb_dotProduct_v_gvec (m : SpikedModel μ n d) (hG : m.GaussianNoise)
    (hd : Tendsto d atTop atTop) :
    TendstoInProb μ (fun N ω => WithLp.ofLp (m.v N) ⬝ᵥ m.gvec N ω) 0 := by
  classical
  have hpsucc : ∀ N, n N = (n N - 1) + 1 := fun N =>
    (Nat.succ_pred_eq_of_pos (m.hn N)).symm
  choose B hW hlaw using fun N => m.exists_block_hasLaw N hG (hpsucc N)
  set r : ∀ N, Fin (d N) → ℝ := fun N a =>
    (Real.sqrt (d N))⁻¹ * WithLp.ofLp (m.v N) a with hrdef
  have hsum : ∀ (N : ℕ) (x : Fin (d N) → ℝ),
      WithLp.ofLp (m.v N) ⬝ᵥ R2.gOf x = ∑ a, r N a * x a := by
    intro N x
    refine Finset.sum_congr rfl fun a _ => ?_
    change WithLp.ofLp (m.v N) a * ((Real.sqrt (d N))⁻¹ * x a) = r N a * x a
    rw [hrdef]
    ring
  have hrsq : ∀ N : ℕ, ∑ a, r N a ^ 2 = ((d N : ℝ))⁻¹ := by
    intro N
    have hterm : ∀ a : Fin (d N),
        r N a ^ 2 = ((d N : ℝ))⁻¹ * (WithLp.ofLp (m.v N) a * WithLp.ofLp (m.v N) a) := by
      intro a
      rw [hrdef]
      change ((Real.sqrt (d N))⁻¹ * WithLp.ofLp (m.v N) a) ^ 2 = _
      rw [mul_pow, ← R2.sqrt_inv_mul_sqrt_inv (d N)]
      ring
    rw [Finset.sum_congr rfl fun a _ => hterm a, ← Finset.mul_sum]
    have hvv : ∑ a, WithLp.ofLp (m.v N) a * WithLp.ofLp (m.v N) a = 1 :=
      R2.dotProduct_ofLp_self (m.hv N)
    rw [hvv, mul_one]
  have hbase : TendstoInProb μ
      (fun N ω => |WithLp.ofLp (m.v N) ⬝ᵥ R2.gOf ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N))|) 0 := by
    refine R2.tendstoInProb_of_integral_sq_le (μ := μ) (pN := fun N => n N - 1) (dN := d)
      (fun N ω => ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N), B N ω)) hlaw
      (fun N x _ => |WithLp.ofLp (m.v N) ⬝ᵥ R2.gOf x|) (fun N => ?_)
      (fun _ _ _ => abs_nonneg _) (fun N _ => ?_) (K := 1) ?_ hd
    · have h1 : Measurable fun q : R2.NoiseSpace (n N - 1) (d N) => ∑ a, r N a * q.1 a :=
        Finset.measurable_sum _ fun a _ =>
          ((measurable_pi_apply a).comp measurable_fst).const_mul _
      have h2 : (fun q : R2.NoiseSpace (n N - 1) (d N) =>
          |WithLp.ofLp (m.v N) ⬝ᵥ R2.gOf q.1|)
          = fun q => |∑ a, r N a * q.1 a| := by
        funext q; rw [hsum]
      rw [h2]
      exact h1.abs
    · have h1 := R2.centered_id.integrable_sq_sum (r N)
      have h2 : (fun x : Fin (d N) → ℝ => |WithLp.ofLp (m.v N) ⬝ᵥ R2.gOf x| ^ 2)
          = fun x => (∑ a, r N a * ((fun t : ℝ => t) (x a))) ^ 2 := by
        funext x; rw [sq_abs, hsum]
      rw [h2]
      exact h1
    · refine Filter.Eventually.of_forall fun N _ => ?_
      have h1 := R2.centered_id.integral_sq_sum (r N)
      have h2 : (fun x : Fin (d N) → ℝ => |WithLp.ofLp (m.v N) ⬝ᵥ R2.gOf x| ^ 2)
          = fun x => (∑ a, r N a * ((fun t : ℝ => t) (x a))) ^ 2 := by
        funext x; rw [sq_abs, hsum]
      rw [h2, h1, R2.integral_sq_gauss, one_mul, hrsq N, one_div]
  refine TendstoInProb.of_le
    (g := fun N ω => |WithLp.ofLp (m.v N) ⬝ᵥ R2.gOf ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N))|)
    (fun N => ?_) hbase
  filter_upwards with ω
  rw [sub_zero, m.gvec_eq_gOf N ω]

/-- `q ⬝ q → θ² + 1`. -/
theorem tendstoInProb_dotProduct_qvec (m : SpikedModel μ n d) (hG : m.GaussianNoise)
    (hd : Tendsto d atTop atTop) :
    TendstoInProb μ (fun N ω => m.qvec N ω ⬝ᵥ m.qvec N ω) (m.θ ^ 2 + 1) := by
  have hfun : ∀ (N : ℕ) (ω : Ω N), m.qvec N ω ⬝ᵥ m.qvec N ω
      = m.θ ^ 2 + 2 * m.θ * (WithLp.ofLp (m.v N) ⬝ᵥ m.gvec N ω)
        + m.gvec N ω ⬝ᵥ m.gvec N ω := by
    intro N ω
    change (m.θ • WithLp.ofLp (m.v N) + m.gvec N ω)
        ⬝ᵥ (m.θ • WithLp.ofLp (m.v N) + m.gvec N ω) = _
    simp only [add_dotProduct, dotProduct_add, smul_dotProduct, dotProduct_smul, smul_eq_mul,
      R2.dotProduct_ofLp_self (m.hv N),
      dotProduct_comm (m.gvec N ω) (WithLp.ofLp (m.v N))]
    ring
  have h := ((TendstoInProb.const μ (m.θ ^ 2)).add
    ((m.tendstoInProb_dotProduct_v_gvec hG hd).const_mul (2 * m.θ))).add
      (m.tendstoInProb_dotProduct_gvec hG hd)
  have hlim : m.θ ^ 2 + 2 * m.θ * 0 + 1 = m.θ ^ 2 + 1 := by ring
  rw [hlim] at h
  simpa only [hfun] using h

variable [∀ N, IsProbabilityMeasure (μ N)] {m : SpikedModel μ n d} {c : ℝ}

set_option linter.unusedSectionVars false

/-- (H2) at the vector `q = θ v + g`: the squared resolvent form converges to
`(θ² + 1) m'(z)` at every real `z` above the bulk edge. -/
theorem tendstoInProb_qform2_qvec (H : m.ResolventLimits c) {z : ℝ} (hz : bulkEdge c < z) :
    TendstoInProb μ (fun N ω => R4.qform2 (m.W0 N ω) z (m.qvec N ω))
      ((m.θ ^ 2 + 1) * MP.mDeriv c z) := by
  have hfun : ∀ (N : ℕ) (ω : Ω N), R4.qform2 (m.W0 N ω) z (m.qvec N ω)
      = m.θ ^ 2 * R4.qform2 (m.W0 N ω) z (WithLp.ofLp (m.v N))
        + 2 * m.θ * R4.cform2 (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω)
        + R4.qform2 (m.W0 N ω) z (m.gvec N ω) := by
    intro N ω
    change R4.qform2 (m.W0 N ω) z (m.θ • WithLp.ofLp (m.v N) + m.gvec N ω) = _
    rw [R4.qform2_smul_add (m.isHermitian_W0 N ω)]
  have h := (((H.vv2 z hz).const_mul (m.θ ^ 2)).add
    ((H.vg2 z hz).const_mul (2 * m.θ))).add (H.gg2 z hz)
  have hlim : m.θ ^ 2 * MP.mDeriv c z + 2 * m.θ * 0 + MP.mDeriv c z
      = (m.θ ^ 2 + 1) * MP.mDeriv c z := by ring
  rw [hlim] at h
  simpa only [hfun] using h

/-! #### The align field, subcritical branch -/

/-- **R6'.** For `θ⁴ ≤ c` the overlap with `v` vanishes in probability, which is the `align`
field of `SingleTableLaw` because `betaSq θ c = 0` there. The real point `z₀` above the bulk
edge is chosen inside the proof: at a fixed `z₀` the bound `8/((θ²+1)² m'(z₀))` need not be
below the target, and `MP.mDeriv_tendsto_atTop` makes it so as `z₀` decreases to the edge. -/
theorem align_tendstoInProb_of_subcritical (H : m.ResolventLimits c) (hc : 0 < c)
    (hreg : m.Regime c) (hG : m.GaussianNoise) (hθ : m.θ ^ 4 ≤ c) :
    TendstoInProb μ (fun N ω => overlap (m.X N ω) (m.v N)) (betaSq m.θ c) := by
  have hbeta : betaSq m.θ c = 0 := by rw [betaSq, if_neg (not_lt.mpr hθ)]
  rw [hbeta]
  intro ε hε
  have hK : (0 : ℝ) < m.θ ^ 2 + 1 := by positivity
  obtain ⟨z₀, hz₀b, -, hz₀K⟩ :=
    R3minus.exists_mDeriv_gt hc (16 / ((m.θ ^ 2 + 1) ^ 2 * ε))
  have hmD : 0 < MP.mDeriv c z₀ := MP.mDeriv_pos hc hz₀b
  have hbound : 8 < ε / 2 * ((m.θ ^ 2 + 1) ^ 2 * MP.mDeriv c z₀) := by
    rw [div_lt_iff₀ (by positivity)] at hz₀K
    nlinarith
  -- the five bad families
  have hB1 : Tendsto (fun N => μ N {ω | gramLamMax (m.X N ω) ≤ z₀}ᶜ) atTop (𝓝 0) := by
    have h := tendsto_measure_lamMax_le_of_subcritical H hc hθ (z₀ - bulkEdge c)
      (by linarith)
    have he : bulkEdge c + (z₀ - bulkEdge c) = z₀ := by ring
    rw [he] at h
    exact tendsto_measure_compl_zero
      (fun N => (m.measurableSet_gramLamMax_le N _).nullMeasurableSet) h
  have hB2 := m.tendstoInProb_dotProduct_qvec hG hreg.2.1 ((m.θ ^ 2 + 1) / 2) (by linarith)
  have hB3 := tendstoInProb_qform2_qvec H hz₀b
    ((m.θ ^ 2 + 1) * MP.mDeriv c z₀ / 2) (by positivity)
  have hB4 := m.tendstoInProb_overlap_rdir hG hreg.2.1 (ε / 4) (by linarith)
  have hB5 : Tendsto (fun N => μ N {ω : Ω N | ¬ TopSimple ((m.X N ω)ᵀ * m.X N ω)
      (isHermitian_transpose_mul_self (m.X N ω))}) atTop (𝓝 0) := by
    have hz : ∀ N, μ N {ω : Ω N | ¬ TopSimple ((m.X N ω)ᵀ * m.X N ω)
        (isHermitian_transpose_mul_self (m.X N ω))} = 0 := fun N =>
      ae_iff.mp (singleTableLaw_topSimple_of_gaussian m hG N)
    simp only [hz]
    exact tendsto_const_nhds
  refine tendsto_measure_zero_of_subset ?_ (tendsto_measure_zero_union hB1
    (tendsto_measure_zero_union hB2 (tendsto_measure_zero_union hB3
      (tendsto_measure_zero_union hB4 hB5))))
  intro N ω hω
  by_contra hcon
  simp only [Set.mem_union, not_or] at hcon
  obtain ⟨g1, g2, g3, g4, g5⟩ := hcon
  have hlam : gramLamMax (m.X N ω) ≤ z₀ := by
    by_contra hx
    exact g1 hx
  have hqq : (m.θ ^ 2 + 1) / 2 < m.qvec N ω ⬝ᵥ m.qvec N ω := by
    have h : ¬ ((m.θ ^ 2 + 1) / 2
        ≤ |m.qvec N ω ⬝ᵥ m.qvec N ω - (m.θ ^ 2 + 1)|) := g2
    have h2 := abs_lt.mp (not_le.mp h)
    linarith [h2.1]
  have hQ2 : (m.θ ^ 2 + 1) * MP.mDeriv c z₀ / 2
      < R4.qform2 (m.W0 N ω) z₀ (m.qvec N ω) := by
    have h : ¬ ((m.θ ^ 2 + 1) * MP.mDeriv c z₀ / 2
        ≤ |R4.qform2 (m.W0 N ω) z₀ (m.qvec N ω) - (m.θ ^ 2 + 1) * MP.mDeriv c z₀|) := g3
    have h2 := abs_lt.mp (not_le.mp h)
    linarith [h2.1]
  have hovr : overlap (m.X N ω) (WithLp.toLp 2 (m.rdir N ω)) < ε / 4 := by
    have h : ¬ (ε / 4
        ≤ |overlap (m.X N ω) (WithLp.toLp 2 (m.rdir N ω)) - 0|) := g4
    have h2 := not_le.mp h
    rw [sub_zero, abs_of_nonneg (overlap_nonneg _ _)] at h2
    exact h2
  have hsimple : TopSimple ((m.X N ω)ᵀ * m.X N ω)
      (isHermitian_transpose_mul_self (m.X N ω)) := not_not.mp g5
  -- the two factors
  have hQ2pos : 0 < R4.qform2 (m.W0 N ω) z₀ (m.qvec N ω) := by
    have : (0 : ℝ) < (m.θ ^ 2 + 1) * MP.mDeriv c z₀ / 2 := by positivity
    linarith
  have hqqpos : 0 < m.qvec N ω ⬝ᵥ m.qvec N ω := by linarith
  have hq := R6.overlap_q_le (m.X N ω) (m.isHermitian_W0 N ω) (m.gram_eq N ω) (m.hd N)
    hsimple hlam hQ2pos
  have hovqnn : 0 ≤ overlap (m.X N ω) (WithLp.toLp 2 (m.qvec N ω)) := overlap_nonneg _ _
  have hovq : overlap (m.X N ω) (WithLp.toLp 2 (m.qvec N ω))
      * R4.qform2 (m.W0 N ω) z₀ (m.qvec N ω) ≤ 1 := by
    rw [le_div_iff₀ hQ2pos] at hq
    linarith [hq]
  have ha2nn : 0 ≤ R6.aOf (WithLp.ofLp (m.v N)) (m.qvec N ω) ^ 2 := sq_nonneg _
  have ha2 : R6.aOf (WithLp.ofLp (m.v N)) (m.qvec N ω) ^ 2
      * (m.qvec N ω ⬝ᵥ m.qvec N ω) ≤ 1 := by
    have h := R6.aOf_sq_mul_le (WithLp.ofLp (m.v N)) (m.qvec N ω)
    rw [R2.dotProduct_ofLp_self (m.hv N)] at h
    exact h
  -- the decomposition
  have hdec := R6.overlap_le_decomp (m.X N ω) (WithLp.ofLp (m.v N)) (m.qvec N ω)
  have hvv : (WithLp.toLp 2 (WithLp.ofLp (m.v N)) : EuclideanSpace ℝ (Fin (d N))) = m.v N := rfl
  rw [hvv] at hdec
  -- the arithmetic
  have hs1 : R6.aOf (WithLp.ofLp (m.v N)) (m.qvec N ω) ^ 2 * (m.θ ^ 2 + 1) ≤ 2 := by
    nlinarith
  have hs2 : overlap (m.X N ω) (WithLp.toLp 2 (m.qvec N ω))
      * ((m.θ ^ 2 + 1) * MP.mDeriv c z₀) ≤ 2 := by
    nlinarith
  have hs3 : (R6.aOf (WithLp.ofLp (m.v N)) (m.qvec N ω) ^ 2
      * overlap (m.X N ω) (WithLp.toLp 2 (m.qvec N ω)))
      * ((m.θ ^ 2 + 1) ^ 2 * MP.mDeriv c z₀) ≤ 4 := by
    have he : (R6.aOf (WithLp.ofLp (m.v N)) (m.qvec N ω) ^ 2
        * overlap (m.X N ω) (WithLp.toLp 2 (m.qvec N ω)))
        * ((m.θ ^ 2 + 1) ^ 2 * MP.mDeriv c z₀)
        = (R6.aOf (WithLp.ofLp (m.v N)) (m.qvec N ω) ^ 2 * (m.θ ^ 2 + 1))
          * (overlap (m.X N ω) (WithLp.toLp 2 (m.qvec N ω))
            * ((m.θ ^ 2 + 1) * MP.mDeriv c z₀)) := by ring
    rw [he]
    nlinarith [mul_nonneg ha2nn hK.le, mul_nonneg hovqnn (mul_pos hK hmD).le]
  have hs4 : 2 * (R6.aOf (WithLp.ofLp (m.v N)) (m.qvec N ω) ^ 2
      * overlap (m.X N ω) (WithLp.toLp 2 (m.qvec N ω))) < ε / 2 := by
    have hpos : (0 : ℝ) < (m.θ ^ 2 + 1) ^ 2 * MP.mDeriv c z₀ := by positivity
    nlinarith
  have hrd : (WithLp.toLp 2 (R6.rvecOf (WithLp.ofLp (m.v N)) (m.qvec N ω))
      : EuclideanSpace ℝ (Fin (d N))) = WithLp.toLp 2 (m.rdir N ω) := rfl
  rw [hrd] at hdec
  have hfin : overlap (m.X N ω) (m.v N) < ε := by linarith
  have hω' : ε ≤ |overlap (m.X N ω) (m.v N) - 0| := hω
  rw [sub_zero, abs_of_nonneg (overlap_nonneg _ _)] at hω'
  linarith

end SpikedModel

end StackedSVD
