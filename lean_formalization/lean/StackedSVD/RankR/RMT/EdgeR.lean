/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.Outliers
import StackedSVD.LinAlg.SpecIdxMeas

/-!
# Task U7, gap G1: the subcritical edge block

Gap G1 of `notes/archive/rankr_plan_A.md` section 4. The mixed case splits the stacked Gram matrix
as `S = W₁ + Q_sup Q_supᵀ` with `W₁ = W₀ + Q_sub Q_subᵀ`, and the count `= s` that
`Frame.specProj_frame_approx` needs comes from `Frame.eigenvalues₀_le_of_split` applied to that
split. The input it asks for is the edge bound for `W₁`, which this file proves.

The route is the plan's, confirmed by the adversarial audit
(`notes/archive/audit_rankr_plan_A_2026-09-01.md`, the verdict rows on G1 and table C2, seed
`2026090101`): `P[lamMax (W₀ + Q_sub Q_subᵀ) ≤ b + ε]` rises to 1.000 with `d` in every cell,
including the tied, the rank-deficient and the exactly critical spikes. G1 is therefore a
convergence in probability, not an almost-sure statement.

Contents:

1. **Deterministic layer** (`EdgeR`). `posSemidef_of_close_to_diag`: a symmetric matrix whose
   entries sit within `δ` of a diagonal matrix with entries at least `δ * t` is positive
   semidefinite (an elementary pairwise bound, `2 |u| |v| ≤ u² + v²`, so no Cauchy-Schwarz).
   `lamMax_le_of_posSemidef`: `z • 1 - A ⪰ 0` bounds the top eigenvalue. `posDef_smul_one_sub`:
   `z • 1 - W₀ ≻ 0` above `lamMax W₀`. `lamMax_add_le_of_posSemidef`: the Schur step, through
   Mathlib's `Matrix.PosDef.fromBlocks₁₁` and `fromBlocks₂₂` at `A = z • 1 - W₀`, `B = Q_sub`,
   `D = 1`. Mind the sign: `resolv W z = (W - z • 1)⁻¹` (`RMT/R4.lean`), so
   `(z • 1 - W₀)⁻¹ = -resolv W₀ z` and the Schur complement `1 - Q_subᵀ (z • 1 - W₀)⁻¹ Q_sub`
   is the plan's `1 + Q_subᵀ G₀(z) Q_sub`.
2. **The transfer** (`RankRStack.tendsto_measure_lamMax_w1R_le`). At `z = bulkEdge c + ε` the
   entries of `1 + Q_subᵀ G₀(z) Q_sub` tend in probability to `diag (1 + (λ_j + 1) m(z))`
   (`ResolventLimitsR`, through `tendstoInProb_cform_qmatR`), and each diagonal limit is at
   least `1 + (√c + 1) m(z) > 0` by `MP.secular_pos_of_subcritical`. That uniform bound is what
   makes the accuracy `δ` a single number for all `t` indices.
3. **The corollary** (`RankRStack.tendsto_measure_eigenvalues₀_le`). With a partition of the
   `r` spikes into `Fin t` subcritical and `Fin u` supercritical ones,
   `Frame.eigenvalues₀_le_of_split` gives `eigenvalues₀ k ≤ bulkEdge c + ε` for every `k ≥ u`,
   with probability tending to 1. This is the hypothesis shape of
   `Frame.norm_sq_specProjTop_split` and `Frame.specProjTop_eq_specProj_add_edge`.

**Zero spikes.** The plan expected a separate branch for `λ_j = 0`, because
`MP.secular_pos_of_subcritical` was read as asking `0 < θ`. It does not: its second argument is
`_hθ : 0 ≤ θ` and it is unused (`RMT/MP.lean:506`). The uniform bound of step 2 replaces the
per-index call altogether, so a zero spike needs no branch. The subcritical hypothesis is
`coreEig j ^ 2 ≤ c`, the same normalization as the supercritical `c < coreEig j ^ 2` of
`RankR/RMT/Outliers.lean`, that is `θ⁴` against `c` with `λ = θ²`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace EdgeR

/-! ### 1. Real symmetry and the two Hermitian bridges -/

section Herm

variable {p q : ℕ}

/-- Over `ℝ`, `Aᴴ = Aᵀ`, so `IsHermitian` reads as symmetry of the transpose. -/
theorem isHermitian_of_transpose_eq {A : Matrix (Fin p) (Fin p) ℝ} (h : Aᵀ = A) :
    A.IsHermitian := by
  change Aᴴ = A
  rw [Matrix.conjTranspose_eq_transpose_of_trivial]
  exact h

/-- The converse bridge. -/
theorem transpose_eq_of_isHermitian {A : Matrix (Fin p) (Fin p) ℝ} (h : A.IsHermitian) :
    Aᵀ = A := by
  rwa [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial] at h

/-- `A Aᵀ` is symmetric. -/
theorem isHermitian_mul_transpose (A : Matrix (Fin p) (Fin q) ℝ) : (A * Aᵀ).IsHermitian := by
  refine isHermitian_of_transpose_eq ?_
  rw [Matrix.transpose_mul, Matrix.transpose_transpose]

end Herm

/-! ### 2. Positive semidefiniteness from closeness to a diagonal matrix

The transfer of step 2 needs a quantitative form. Write `M = D + E` with `D = diag L` and
`|E a b| ≤ δ`. Then `|x ⬝ᵥ (E x)| ≤ δ t ∑ x²` by `2 |u| |v| ≤ u² + v²`, while
`x ⬝ᵥ (D x) ≥ (min L) ∑ x²`. No Cauchy-Schwarz is needed. -/

section CloseToDiag

variable {t : ℕ}

/-- **The PSD transfer.** A symmetric matrix within `δ` of `diag L`, entrywise, with every
`L a` at least `δ * t`, is positive semidefinite. -/
theorem posSemidef_of_close_to_diag {M : Matrix (Fin t) (Fin t) ℝ} (hsymm : Mᵀ = M)
    {L : Fin t → ℝ} {δ : ℝ} (hδ : 0 ≤ δ) (hL : ∀ a, δ * t ≤ L a)
    (hclose : ∀ a b, |M a b - (if a = b then L a else 0)| ≤ δ) :
    M.PosSemidef := by
  refine Matrix.PosSemidef.of_dotProduct_mulVec_nonneg (isHermitian_of_transpose_eq hsymm)
    fun x => ?_
  have hstar : star x = x := rfl
  rw [hstar, Matrix.dot_mulVec_eq_sum_sum]
  set y : ℝ := ∑ a : Fin t, x a ^ 2 with hy
  have hy0 : 0 ≤ y := Finset.sum_nonneg fun a _ => sq_nonneg _
  -- the entrywise error bound, symmetric in the two indices
  have herrbound : ∀ b a : Fin t, |x a * (M a b - (if a = b then L a else 0)) * x b|
      ≤ δ * ((x a ^ 2 + x b ^ 2) / 2) := by
    intro b a
    have h1 : |x a * (M a b - (if a = b then L a else 0)) * x b|
        = |x a| * |M a b - (if a = b then L a else 0)| * |x b| := by
      rw [abs_mul, abs_mul]
    rw [h1]
    have h3 : |x a| * |M a b - (if a = b then L a else 0)| * |x b| ≤ |x a| * δ * |x b| := by
      have h2 : |x a| * |M a b - (if a = b then L a else 0)| ≤ |x a| * δ :=
        mul_le_mul_of_nonneg_left (hclose a b) (abs_nonneg _)
      exact mul_le_mul_of_nonneg_right h2 (abs_nonneg _)
    refine h3.trans ?_
    have h4 : 2 * |x a| * |x b| ≤ |x a| ^ 2 + |x b| ^ 2 := two_mul_le_add_sq _ _
    rw [sq_abs, sq_abs] at h4
    have h5 : δ * (2 * |x a| * |x b|) ≤ δ * (x a ^ 2 + x b ^ 2) :=
      mul_le_mul_of_nonneg_left h4 hδ
    nlinarith [h5]
  -- column by column
  have hstep : ∀ b : Fin t, L b * x b ^ 2 - δ / 2 * (y + (t : ℝ) * x b ^ 2)
      ≤ ∑ a, x a * M a b * x b := by
    intro b
    have hb1 : ∑ a, x a * M a b * x b
        = (∑ a, x a * (if a = b then L a else 0) * x b)
          + ∑ a, x a * (M a b - (if a = b then L a else 0)) * x b := by
      rw [← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl fun a _ => by ring
    have hb2 : (∑ a, x a * (if a = b then L a else 0) * x b) = L b * x b ^ 2 := by
      rw [Finset.sum_eq_single b]
      · rw [if_pos rfl]; ring
      · intro a _ hab
        rw [if_neg hab]; ring
      · intro hb
        exact absurd (Finset.mem_univ b) hb
    have hb3 : |∑ a, x a * (M a b - (if a = b then L a else 0)) * x b|
        ≤ ∑ a : Fin t, δ * ((x a ^ 2 + x b ^ 2) / 2) :=
      (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun a _ => herrbound b a)
    have hb4 : ∑ a : Fin t, δ * ((x a ^ 2 + x b ^ 2) / 2)
        = δ / 2 * (y + (t : ℝ) * x b ^ 2) := by
      have hsp : ∑ a : Fin t, δ * ((x a ^ 2 + x b ^ 2) / 2)
          = (∑ a : Fin t, δ / 2 * x a ^ 2) + ∑ _a : Fin t, δ / 2 * x b ^ 2 := by
        rw [← Finset.sum_add_distrib]
        exact Finset.sum_congr rfl fun a _ => by ring
      rw [hsp, ← Finset.mul_sum, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
        nsmul_eq_mul, ← hy]
      ring
    rw [hb1, hb2]
    have hlow := (abs_le.mp (hb3.trans_eq hb4)).1
    linarith
  have h6 : ∑ b, (L b * x b ^ 2 - δ / 2 * (y + (t : ℝ) * x b ^ 2))
      ≤ ∑ b, ∑ a, x a * M a b * x b :=
    Finset.sum_le_sum fun b _ => hstep b
  have h7 : δ * (t : ℝ) * y ≤ ∑ b, L b * x b ^ 2 := by
    have hrw : δ * (t : ℝ) * y = ∑ b : Fin t, δ * (t : ℝ) * x b ^ 2 := by
      rw [hy, Finset.mul_sum]
    rw [hrw]
    exact Finset.sum_le_sum fun b _ => mul_le_mul_of_nonneg_right (hL b) (sq_nonneg _)
  have h8 : ∑ b, (L b * x b ^ 2 - δ / 2 * (y + (t : ℝ) * x b ^ 2))
      = (∑ b, L b * x b ^ 2) - δ * (t : ℝ) * y := by
    rw [Finset.sum_sub_distrib]
    congr 1
    have hsp : ∑ b : Fin t, δ / 2 * (y + (t : ℝ) * x b ^ 2)
        = (∑ _b : Fin t, δ / 2 * y) + ∑ b : Fin t, δ / 2 * (t : ℝ) * x b ^ 2 := by
      rw [← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl fun b _ => by ring
    rw [hsp, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul,
      ← Finset.mul_sum, ← hy]
    ring
  linarith

end CloseToDiag

/-! ### 3. The top eigenvalue and the Loewner order -/

section LamMax

variable {p : ℕ} {A W₀ : Matrix (Fin p) (Fin p) ℝ} {z : ℝ}

/-- The Rayleigh bound: `x ⬝ᵥ (A x) ≤ lamMax A * (x ⬝ᵥ x)`. -/
theorem dotProduct_le_lamMax (hA : A.IsHermitian) (x : Fin p → ℝ) :
    x ⬝ᵥ (A *ᵥ x) ≤ lamMax A hA * (x ⬝ᵥ x) := by
  have h1 : x ⬝ᵥ (A *ᵥ x) = ∑ a, hA.eigenvalues a * ((R4.eigU hA)ᵀ *ᵥ x) a ^ 2 := by
    conv_lhs => rw [← R4.eigU_conj hA]
    exact R4.dotProduct_conj _ _ _
  have h2 : x ⬝ᵥ x = ∑ a, ((R4.eigU hA)ᵀ *ᵥ x) a ^ 2 := by
    rw [← R4.dotProduct_transpose_eigU hA x]
    simp only [dotProduct]
    exact Finset.sum_congr rfl fun a _ => (sq _).symm
  rw [h1, h2, Finset.mul_sum]
  refine Finset.sum_le_sum fun a _ => ?_
  have h3 := R4.eigenvalues_le_lamMax hA a
  nlinarith [sq_nonneg (((R4.eigU hA)ᵀ *ᵥ x) a)]

/-- **`z • 1 - W₀ ≻ 0` above the top eigenvalue.** -/
theorem posDef_smul_one_sub (hW₀ : W₀.IsHermitian) (hz : lamMax W₀ hW₀ < z) :
    (z • (1 : Matrix (Fin p) (Fin p) ℝ) - W₀).PosDef := by
  have hherm : (z • (1 : Matrix (Fin p) (Fin p) ℝ) - W₀).IsHermitian := by
    refine isHermitian_of_transpose_eq ?_
    rw [Matrix.transpose_sub, Matrix.transpose_smul, Matrix.transpose_one,
      transpose_eq_of_isHermitian hW₀]
  refine Matrix.PosDef.of_dotProduct_mulVec_pos hherm fun x hx => ?_
  have hstar : star x = x := rfl
  rw [hstar, Matrix.sub_mulVec, dotProduct_sub, Matrix.smul_mulVec, Matrix.one_mulVec,
    dotProduct_smul, smul_eq_mul]
  have hnn : (0 : ℝ) ≤ x ⬝ᵥ x := Frame.dot_self_nonneg x
  have hxx : 0 < x ⬝ᵥ x := by
    rcases hnn.lt_or_eq with hlt | heq
    · exact hlt
    · exact absurd (dotProduct_self_eq_zero.mp heq.symm) hx
  have h1 := dotProduct_le_lamMax hW₀ x
  nlinarith [h1, hxx]

/-- **From the Loewner order to the top eigenvalue.** -/
theorem lamMax_le_of_posSemidef (hp : 0 < p) (hA : A.IsHermitian)
    (h : (z • (1 : Matrix (Fin p) (Fin p) ℝ) - A).PosSemidef) : lamMax A hA ≤ z := by
  classical
  obtain ⟨i, hi⟩ := R4.exists_eigenvalues_eq_lamMax hA hp
  set e : Fin p → ℝ := fun a => if a = i then 1 else 0 with he
  set x : Fin p → ℝ := R4.eigU hA *ᵥ e with hx
  have hUx : (R4.eigU hA)ᵀ *ᵥ x = e := by
    rw [hx, Matrix.mulVec_mulVec, R4.transpose_eigU_mul, Matrix.one_mulVec]
  have hsq : ∀ a : Fin p, e a ^ 2 = if a = i then 1 else 0 := by
    intro a
    simp only [he]
    by_cases hai : a = i <;> simp [hai]
  have h1 : x ⬝ᵥ (A *ᵥ x) = hA.eigenvalues i := by
    have hexp : x ⬝ᵥ (A *ᵥ x) = ∑ a, hA.eigenvalues a * ((R4.eigU hA)ᵀ *ᵥ x) a ^ 2 := by
      conv_lhs => rw [← R4.eigU_conj hA]
      exact R4.dotProduct_conj _ _ _
    rw [hexp, hUx]
    simp [hsq]
  have h2 : x ⬝ᵥ x = 1 := by
    have hexp : x ⬝ᵥ x = ∑ a, ((R4.eigU hA)ᵀ *ᵥ x) a ^ 2 := by
      rw [← R4.dotProduct_transpose_eigU hA x]
      simp only [dotProduct]
      exact Finset.sum_congr rfl fun a _ => (sq _).symm
    rw [hexp, hUx]
    simp [hsq]
  have hnn := h.dotProduct_mulVec_nonneg x
  have hstar : star x = x := rfl
  rw [hstar, Matrix.sub_mulVec, dotProduct_sub, Matrix.smul_mulVec, Matrix.one_mulVec,
    dotProduct_smul, smul_eq_mul, h1, h2, mul_one] at hnn
  rw [← hi]
  linarith

end LamMax

/-! ### 4. The Schur step -/

section Schur

variable {p t : ℕ}

/-- The `(a, b)` entry of `Qᵀ G₀(z) Q` is the bilinear form of the two columns. -/
theorem cform_eq_entry (W : Matrix (Fin p) (Fin p) ℝ) (z : ℝ)
    (Qs : Matrix (Fin p) (Fin t) ℝ) (a b : Fin t) :
    (Qsᵀ * R4.resolv W z * Qs) a b
      = R4.cform W z (fun k => Qs k a) (fun k => Qs k b) := by
  have hL : (Qsᵀ * R4.resolv W z * Qs) a b
      = ∑ l, ∑ k, Qs k a * R4.resolv W z k l * Qs l b := by
    rw [Matrix.mul_apply]
    refine Finset.sum_congr rfl fun l _ => ?_
    rw [Matrix.mul_apply, Finset.sum_mul]
    exact Finset.sum_congr rfl fun k _ => by rw [Matrix.transpose_apply]
  have hR : R4.cform W z (fun k => Qs k a) (fun k => Qs k b)
      = ∑ k, ∑ l, Qs k a * R4.resolv W z k l * Qs l b := by
    rw [R4.cform]
    simp only [dotProduct]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [Matrix.mulVec_apply_eq_sum, Finset.mul_sum]
    exact Finset.sum_congr rfl fun l _ => by ring
  rw [hL, hR, Finset.sum_comm]

/-- **The Schur step.** Above `lamMax W₀`, positive semidefiniteness of
`1 + Q_subᵀ G₀(z) Q_sub` bounds the top eigenvalue of `W₀ + Q_sub Q_subᵀ` by `z`.

Mathlib's `Matrix.PosDef.fromBlocks₁₁` and `fromBlocks₂₂` are the two halves, at
`A = z • 1 - W₀`, `B = Q_sub`, `D = 1`. The sign enters through
`(z • 1 - W₀)⁻¹ = -resolv W₀ z`. -/
theorem lamMax_add_le_of_posSemidef (hp : 0 < p) {W₀ : Matrix (Fin p) (Fin p) ℝ}
    (hW₀ : W₀.IsHermitian) {Qs : Matrix (Fin p) (Fin t) ℝ} {z : ℝ}
    (hz : lamMax W₀ hW₀ < z)
    (hpsd : ((1 : Matrix (Fin t) (Fin t) ℝ) + Qsᵀ * R4.resolv W₀ z * Qs).PosSemidef)
    (hW₁ : (W₀ + Qs * Qsᵀ).IsHermitian) :
    lamMax (W₀ + Qs * Qsᵀ) hW₁ ≤ z := by
  classical
  have hA : (z • (1 : Matrix (Fin p) (Fin p) ℝ) - W₀).PosDef := posDef_smul_one_sub hW₀ hz
  have hinvA : Invertible (z • (1 : Matrix (Fin p) (Fin p) ℝ) - W₀) := hA.isUnit.invertible
  have hinv1 : Invertible (1 : Matrix (Fin t) (Fin t) ℝ) := ⟨1, one_mul 1, mul_one 1⟩
  have hinv : (z • (1 : Matrix (Fin p) (Fin p) ℝ) - W₀)⁻¹ = -R4.resolv W₀ z := by
    refine Matrix.inv_eq_right_inv ?_
    have h1 : z • (1 : Matrix (Fin p) (Fin p) ℝ) - W₀ = -(W₀ - z • 1) := by
      rw [neg_sub]
    rw [h1, neg_mul_neg]
    exact R4.mul_resolv hW₀ hz
  have hCT : (Qs)ᴴ = (Qs)ᵀ := by
    rw [Matrix.conjTranspose_eq_transpose_of_trivial]
  have key : ((1 : Matrix (Fin t) (Fin t) ℝ)
      - Qsᴴ * (z • (1 : Matrix (Fin p) (Fin p) ℝ) - W₀)⁻¹ * Qs).PosSemidef := by
    rw [hCT, hinv]
    have hrw : (1 : Matrix (Fin t) (Fin t) ℝ) - Qsᵀ * (-R4.resolv W₀ z) * Qs
        = (1 : Matrix (Fin t) (Fin t) ℝ) + Qsᵀ * R4.resolv W₀ z * Qs := by
      rw [Matrix.mul_neg, Matrix.neg_mul, sub_neg_eq_add]
    rw [hrw]
    exact hpsd
  have h₁ := Matrix.PosDef.fromBlocks₁₁ (A := z • (1 : Matrix (Fin p) (Fin p) ℝ) - W₀) Qs
    (1 : Matrix (Fin t) (Fin t) ℝ) hA
  have h₂ := Matrix.PosDef.fromBlocks₂₂ (z • (1 : Matrix (Fin p) (Fin p) ℝ) - W₀) Qs
    (D := (1 : Matrix (Fin t) (Fin t) ℝ)) Matrix.PosDef.one
  have hres := h₂.mp (h₁.mpr key)
  have hgoal : (z • (1 : Matrix (Fin p) (Fin p) ℝ) - (W₀ + Qs * Qsᵀ)).PosSemidef := by
    have hrw : z • (1 : Matrix (Fin p) (Fin p) ℝ) - W₀
          - Qs * (1 : Matrix (Fin t) (Fin t) ℝ)⁻¹ * Qsᴴ
        = z • (1 : Matrix (Fin p) (Fin p) ℝ) - (W₀ + Qs * Qsᵀ) := by
      rw [inv_one, Matrix.mul_one, hCT, sub_sub]
    rwa [hrw] at hres
  exact lamMax_le_of_posSemidef hp hW₁ hgoal

/-- Splitting `Q Qᵀ` along a partition of the columns into two blocks. -/
theorem mul_transpose_split {rr tt uu : ℕ} (Q : Matrix (Fin p) (Fin rr) ℝ)
    (e : Fin tt ⊕ Fin uu ≃ Fin rr) :
    Q * Qᵀ
      = Q.submatrix id (fun a => e (Sum.inl a)) * (Q.submatrix id (fun a => e (Sum.inl a)))ᵀ
        + Q.submatrix id (fun b => e (Sum.inr b))
          * (Q.submatrix id (fun b => e (Sum.inr b)))ᵀ := by
  ext k l
  simp only [Matrix.add_apply, Matrix.mul_apply, Matrix.transpose_apply,
    Matrix.submatrix_apply, id_eq]
  have h1 : ∑ j : Fin rr, Q k j * Q l j
      = ∑ w : Fin tt ⊕ Fin uu, Q k (e w) * Q l (e w) :=
    (Equiv.sum_comp e fun j => Q k j * Q l j).symm
  rw [h1, Fintype.sum_sum_type]

end Schur

end EdgeR

/-! ### 5. The model layer: the subcritical block `W₁` and its edge -/

namespace RankRStack

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-! ### 5a. Measurability of the table and of the sorted eigenvalues

`RankR/RMT/AlignOutG.lean` and `RankR/RMT/EdgeGlueR.lean` each carried a copy of
`measurable_X` and of `measurableSet_eigenvalues₀_le`; the two files are siblings, so the
second cleanup pass (2026-09-02) moved one copy of each here, the ancestor both import. -/

/-- The table of a `RankRStack` is measurable, entrywise: the signal part is constant in `ω`
and the noise is measurable (`measurable_E`, `RankR/RMT/Stack.lean`). The rank-1 mirror is
`SpikedModel.measurable_X` (`RMT/R3minus.lean`). -/
theorem measurable_X (s : RankRStack μ ns d r) (N : ℕ) : Measurable (s.X N) := by
  refine measurable_pi_lambda _ fun q => measurable_pi_lambda _ fun l => ?_
  have h : (fun ω => s.X N ω q l)
      = fun ω => s.signalPart N q l + s.E N ω q l := by
    funext ω
    rw [s.X_eq N ω]
    rfl
  rw [h]
  exact measurable_const.add
    ((measurable_pi_apply l).comp ((measurable_pi_apply q).comp (s.measurable_E N)))

/-- The sorted eigenvalue of the Gram matrix at a fixed index, as a function of `ω`. -/
theorem measurable_gram_eigenvalues₀ (s : RankRStack μ ns d r) (N : ℕ)
    (k : Fin (Fintype.card (Fin (d N)))) :
    Measurable fun ω => (s.isHermitian_gram N ω).eigenvalues₀ k := by
  have h : (fun ω => (s.isHermitian_gram N ω).eigenvalues₀ k)
      = fun ω => gramEig (s.X N ω) (k : ℕ) := by
    funext ω
    rw [gramEig_of_lt (s.X N ω) k.isLt]
    rfl
  rw [h]
  exact (measurable_gramEig (k : ℕ)).comp (s.measurable_X N)

/-- The count event of gap G1 (`tendsto_measure_eigenvalues₀_le` below) is measurable, so
`tendsto_measure_compl_zero` applies to it. It is a finite intersection of level sets of the
sorted eigenvalues. -/
theorem measurableSet_eigenvalues₀_le (s : RankRStack μ ns d r) (N : ℕ) (u : ℕ) (x : ℝ) :
    MeasurableSet {ω : Ω N | ∀ k : Fin (Fintype.card (Fin (d N))), u ≤ (k : ℕ) →
      (s.isHermitian_gram N ω).eigenvalues₀ k ≤ x} := by
  classical
  have hset : {ω : Ω N | ∀ k : Fin (Fintype.card (Fin (d N))), u ≤ (k : ℕ) →
        (s.isHermitian_gram N ω).eigenvalues₀ k ≤ x}
      = ⋂ k : Fin (Fintype.card (Fin (d N))),
          {ω : Ω N | u ≤ (k : ℕ) → (s.isHermitian_gram N ω).eigenvalues₀ k ≤ x} := by
    ext ω
    simp only [Set.mem_ofPred_eq, Set.mem_iInter]
  rw [hset]
  refine MeasurableSet.iInter fun k => ?_
  by_cases hk : u ≤ (k : ℕ)
  · have he : {ω : Ω N | u ≤ (k : ℕ) → (s.isHermitian_gram N ω).eigenvalues₀ k ≤ x}
        = (fun ω => (s.isHermitian_gram N ω).eigenvalues₀ k) ⁻¹' Set.Iic x := by
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_preimage, Set.mem_Iic]
      exact ⟨fun h => h hk, fun h _ => h⟩
    rw [he]
    exact s.measurable_gram_eigenvalues₀ N k measurableSet_Iic
  · have he : {ω : Ω N | u ≤ (k : ℕ) → (s.isHermitian_gram N ω).eigenvalues₀ k ≤ x}
        = Set.univ := by
      ext ω
      simp only [Set.mem_ofPred_eq, Set.mem_univ, iff_true]
      exact fun hcon => absurd hcon hk
    rw [he]
    exact MeasurableSet.univ

/-- `Q_sub`, the column submatrix of `Q` at the index map `f`. -/
noncomputable def qsubR (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (ns N)) (Fin r) ℝ) {t : ℕ} (f : Fin t → Fin r) :
    Matrix (Fin (d N)) (Fin t) ℝ :=
  (s.qmatR N ω U).submatrix id f

theorem qsubR_apply (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (ns N)) (Fin r) ℝ) {t : ℕ} (f : Fin t → Fin r)
    (k : Fin (d N)) (a : Fin t) :
    s.qsubR N ω U f k a = s.qmatR N ω U k (f a) := rfl

/-- `W₁ = W₀ + Q_sub Q_subᵀ`, the block whose edge closes the count of gap G1. -/
noncomputable def w1R (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (ns N)) (Fin r) ℝ) {t : ℕ} (f : Fin t → Fin r) :
    Matrix (Fin (d N)) (Fin (d N)) ℝ :=
  s.rankRW0 N ω U + s.qsubR N ω U f * (s.qsubR N ω U f)ᵀ

theorem isHermitian_w1R (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (ns N)) (Fin r) ℝ) {t : ℕ} (f : Fin t → Fin r) :
    (s.w1R N ω U f).IsHermitian :=
  (s.isHermitian_rankRW0 N ω U).add (EdgeR.isHermitian_mul_transpose _)

/-- **Gap G1** (`notes/archive/rankr_plan_A.md` section 4). With every column of `Q_sub`
subcritical, `lamMax (W₀ + Q_sub Q_subᵀ) ≤ bulkEdge c + ε` with probability tending to 1.

Numeric validation: `notes/archive/audit_rankr_plan_A_2026-09-01.md` table C2, seed `2026090101`.
The index map `f` selects the subcritical columns; it must be injective, otherwise the limit
matrix is not diagonal. -/
theorem tendsto_measure_lamMax_w1R_le [∀ N, IsProbabilityMeasure (μ N)]
    (s : RankRStack μ ns d r) {U : (N : ℕ) → Matrix (Fin (ns N)) (Fin r) ℝ} {c : ℝ}
    (hc : 0 < c)
    (h : ResolventLimitsR μ (fun N ω => s.rankRW0 N ω (U N))
      (fun N ω => s.isHermitian_rankRW0 N ω (U N)) s.spikeVec
      (fun k N ω => fun l => ((s.E N ω)ᵀ * U N) l k) c)
    {t : ℕ} {f : Fin t → Fin r} (hf : Function.Injective f)
    (hsub : ∀ a : Fin t, s.coreEig (f a) ^ 2 ≤ c) :
    ∀ ε > 0, Tendsto (fun N => μ N
        {ω | lamMax (s.w1R N ω (U N) f) (s.isHermitian_w1R N ω (U N) f)
          ≤ bulkEdge c + ε}) atTop (𝓝 1) := by
  classical
  intro ε hε
  have hbz : bulkEdge c < bulkEdge c + ε := by linarith
  -- the uniform lower bound on the diagonal limit, at the worst subcritical spike
  have hsc0 : (0 : ℝ) ≤ Real.sqrt c := Real.sqrt_nonneg c
  have hL0pos : 0 < 1 + (Real.sqrt c + 1) * MP.m c (bulkEdge c + ε) := by
    have h4 : Real.sqrt (Real.sqrt c) ^ 4 ≤ c := by
      rw [OutliersR.sqrt_pow_four hsc0, Real.sq_sqrt hc.le]
    have hpos := MP.secular_pos_of_subcritical hc (Real.sqrt_nonneg (Real.sqrt c)) h4 hbz
    rwa [Real.sq_sqrt hsc0] at hpos
  have hmneg : MP.m c (bulkEdge c + ε) < 0 := (MP.m_mem_Ioo hc hbz).2
  set L : Fin t → ℝ := fun a => 1 + (s.coreEig (f a) + 1) * MP.m c (bulkEdge c + ε) with hLdef
  have hLle : ∀ a, 1 + (Real.sqrt c + 1) * MP.m c (bulkEdge c + ε) ≤ L a := by
    intro a
    have hle : s.coreEig (f a) ≤ Real.sqrt c := by
      have h1 : Real.sqrt (s.coreEig (f a) ^ 2) ≤ Real.sqrt c := Real.sqrt_le_sqrt (hsub a)
      rwa [Real.sqrt_sq (s.coreEig_nonneg (f a))] at h1
    rw [hLdef]
    nlinarith [hmneg.le, hle]
  -- the single accuracy
  have ht1 : (0 : ℝ) < (t : ℝ) + 1 := by positivity
  set δ : ℝ := (1 + (Real.sqrt c + 1) * MP.m c (bulkEdge c + ε)) / ((t : ℝ) + 1) with hδdef
  have hδpos : 0 < δ := div_pos hL0pos ht1
  have hδt : ∀ a, δ * (t : ℝ) ≤ L a := by
    intro a
    refine le_trans ?_ (hLle a)
    rw [hδdef, div_mul_eq_mul_div, div_le_iff₀ ht1]
    nlinarith [hL0pos.le, Nat.cast_nonneg (α := ℝ) t]
  -- the two bad families
  have hε2 : (0 : ℝ) < ε / 2 := by linarith
  have hedgeC : Tendsto (fun N => μ N
      {ω | lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
        ≤ bulkEdge c + ε / 2}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (s.measurableSet_lamMax_rankRW0_le N (U N)
        (bulkEdge c + ε / 2)).nullMeasurableSet)
      (h.edge (ε / 2) hε2)
  have hET : ∀ q : Fin t × Fin t, Tendsto (fun N => μ N
      {ω | δ ≤ |R4.cform (s.rankRW0 N ω (U N)) (bulkEdge c + ε)
        (fun l => s.qmatR N ω (U N) l (f q.1)) (fun l => s.qmatR N ω (U N) l (f q.2))
        - (if f q.1 = f q.2 then (s.coreEig (f q.2) + 1)
            * MP.m c (bulkEdge c + ε) else 0)|}) atTop (𝓝 0) := fun q =>
    s.tendstoInProb_cform_qmatR h (f q.1) (f q.2) hbz δ hδpos
  have hzero : Tendsto (fun N => μ N
      ({ω | lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
            ≤ bulkEdge c + ε / 2}ᶜ
        ∪ (⋃ q : Fin t × Fin t, {ω | δ ≤ |R4.cform (s.rankRW0 N ω (U N)) (bulkEdge c + ε)
            (fun l => s.qmatR N ω (U N) l (f q.1)) (fun l => s.qmatR N ω (U N) l (f q.2))
            - (if f q.1 = f q.2 then (s.coreEig (f q.2) + 1)
                * MP.m c (bulkEdge c + ε) else 0)|}))) atTop (𝓝 0) :=
    tendsto_measure_zero_union hedgeC (tendsto_measure_zero_iUnion hET)
  have hincl : ∀ N, ({ω | lamMax (s.w1R N ω (U N) f) (s.isHermitian_w1R N ω (U N) f)
        ≤ bulkEdge c + ε} : Set (Ω N))ᶜ
      ⊆ {ω | lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
            ≤ bulkEdge c + ε / 2}ᶜ
        ∪ (⋃ q : Fin t × Fin t, {ω | δ ≤ |R4.cform (s.rankRW0 N ω (U N)) (bulkEdge c + ε)
            (fun l => s.qmatR N ω (U N) l (f q.1)) (fun l => s.qmatR N ω (U N) l (f q.2))
            - (if f q.1 = f q.2 then (s.coreEig (f q.2) + 1)
                * MP.m c (bulkEdge c + ε) else 0)|}) := by
    intro N ω hω
    by_contra hbad
    have hedge : lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
        ≤ bulkEdge c + ε / 2 := by
      by_contra hxx
      exact hbad (Set.mem_union_left _ hxx)
    have hclose : ∀ q : Fin t × Fin t,
        |R4.cform (s.rankRW0 N ω (U N)) (bulkEdge c + ε)
          (fun l => s.qmatR N ω (U N) l (f q.1)) (fun l => s.qmatR N ω (U N) l (f q.2))
          - (if f q.1 = f q.2 then (s.coreEig (f q.2) + 1)
              * MP.m c (bulkEdge c + ε) else 0)| ≤ δ := by
      intro q
      by_contra hxx
      exact hbad (Set.mem_union_right _ (Set.mem_iUnion.mpr ⟨q, (not_le.mp hxx).le⟩))
    have hzlam : lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
        < bulkEdge c + ε := by linarith
    have hM : ((1 : Matrix (Fin t) (Fin t) ℝ)
        + (s.qsubR N ω (U N) f)ᵀ * R4.resolv (s.rankRW0 N ω (U N)) (bulkEdge c + ε)
          * s.qsubR N ω (U N) f).PosSemidef := by
      refine EdgeR.posSemidef_of_close_to_diag ?_ hδpos.le hδt ?_
      · rw [Matrix.transpose_add, Matrix.transpose_one, Matrix.transpose_mul,
          Matrix.transpose_mul, Matrix.transpose_transpose,
          R4.transpose_resolv (s.isHermitian_rankRW0 N ω (U N)), Matrix.mul_assoc]
      · intro a b
        have hentry : ((1 : Matrix (Fin t) (Fin t) ℝ)
              + (s.qsubR N ω (U N) f)ᵀ * R4.resolv (s.rankRW0 N ω (U N)) (bulkEdge c + ε)
                * s.qsubR N ω (U N) f) a b
            - (if a = b then L a else 0)
            = R4.cform (s.rankRW0 N ω (U N)) (bulkEdge c + ε)
                (fun l => s.qmatR N ω (U N) l (f a)) (fun l => s.qmatR N ω (U N) l (f b))
              - (if f a = f b then (s.coreEig (f b) + 1)
                  * MP.m c (bulkEdge c + ε) else 0) := by
          have hcol : ∀ jj : Fin t, (fun k => s.qsubR N ω (U N) f k jj)
              = fun l => s.qmatR N ω (U N) l (f jj) := fun _ => rfl
          rw [Matrix.add_apply, EdgeR.cform_eq_entry, hcol a, hcol b, Matrix.one_apply, hLdef]
          rcases eq_or_ne a b with rfl | hab
          · rw [if_pos rfl, if_pos rfl, if_pos rfl]
            ring
          · rw [if_neg hab, if_neg hab, if_neg (fun hcon => hab (hf hcon))]
            ring
        rw [hentry]
        exact hclose (a, b)
    exact hω (EdgeR.lamMax_add_le_of_posSemidef (s.hd N)
      (s.isHermitian_rankRW0 N ω (U N)) hzlam hM (s.isHermitian_w1R N ω (U N) f))
  exact tendsto_measure_one_of_bad hincl hzero

/-- **The count of gap G1.** With the `r` spikes partitioned into `Fin t` subcritical ones and
`Fin u` supercritical ones, every sorted eigenvalue of the Gram matrix at index `u` or above is
at most `bulkEdge c + ε`, with probability tending to 1.

This is the hypothesis shape of `Frame.norm_sq_specProjTop_split` and of
`Frame.specProjTop_eq_specProj_add_edge`, which task U7c consumes. -/
theorem tendsto_measure_eigenvalues₀_le [∀ N, IsProbabilityMeasure (μ N)]
    (s : RankRStack μ ns d r) {U : (N : ℕ) → Matrix (Fin (ns N)) (Fin r) ℝ} {c : ℝ}
    (hc : 0 < c)
    (h : ResolventLimitsR μ (fun N ω => s.rankRW0 N ω (U N))
      (fun N ω => s.isHermitian_rankRW0 N ω (U N)) s.spikeVec
      (fun k N ω => fun l => ((s.E N ω)ᵀ * U N) l k) c)
    (hgram : ∀ N ω, s.gram N ω
      = s.rankRW0 N ω (U N) + s.qmatR N ω (U N) * (s.qmatR N ω (U N))ᵀ)
    {t u : ℕ} (e : Fin t ⊕ Fin u ≃ Fin r)
    (hsub : ∀ a : Fin t, s.coreEig (e (Sum.inl a)) ^ 2 ≤ c) :
    ∀ ε > 0, Tendsto (fun N => μ N
        {ω | ∀ k : Fin (Fintype.card (Fin (d N))), u ≤ (k : ℕ) →
          (s.isHermitian_gram N ω).eigenvalues₀ k ≤ bulkEdge c + ε}) atTop (𝓝 1) := by
  classical
  intro ε hε
  have hf : Function.Injective (fun a : Fin t => e (Sum.inl a)) :=
    e.injective.comp Sum.inl_injective
  have hmain := s.tendsto_measure_lamMax_w1R_le hc h hf hsub ε hε
  have hsubset : ∀ N, {ω | lamMax (s.w1R N ω (U N) (fun a : Fin t => e (Sum.inl a)))
        (s.isHermitian_w1R N ω (U N) (fun a : Fin t => e (Sum.inl a))) ≤ bulkEdge c + ε}
      ⊆ {ω : Ω N | ∀ k : Fin (Fintype.card (Fin (d N))), u ≤ (k : ℕ) →
          (s.isHermitian_gram N ω).eigenvalues₀ k ≤ bulkEdge c + ε} := by
    intro N ω hω k hk
    have hSeq : s.gram N ω
        = s.w1R N ω (U N) (fun a : Fin t => e (Sum.inl a))
          + s.qsubR N ω (U N) (fun b : Fin u => e (Sum.inr b))
            * (s.qsubR N ω (U N) (fun b : Fin u => e (Sum.inr b)))ᵀ := by
      rw [hgram N ω, EdgeR.mul_transpose_split (s.qmatR N ω (U N)) e, w1R, qsubR, qsubR,
        add_assoc]
    exact Frame.eigenvalues₀_le_of_split
      (s.isHermitian_w1R N ω (U N) (fun a : Fin t => e (Sum.inl a)))
      (s.isHermitian_gram N ω) hSeq hω k hk
  exact tendsto_of_tendsto_of_tendsto_of_le_of_le hmain tendsto_const_nhds
    (fun N => measure_mono (hsubset N)) (fun N => prob_le_one)

end RankRStack

end StackedSVD
