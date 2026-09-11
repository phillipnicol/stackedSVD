/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.Flatten
import StackedSVD.Scalars

/-!
# The block-diagonal Woodbury identity of the exactly aligned model

Task D3det of the D32 campaign, item 1 of `paper_edits.md` section E4. The corollary
`thm:rank_r_svdstack` (`main_paper.tex:2391-2414`) reads the rank-`r` svdstack limit off the
scalars

    β_ij² = (θ_ij⁴ - c_i) / (θ_ij⁴ + θ_ij²) · 1{θ_ij⁴ ≥ c_i},   S_j = Σ_i β_ij² / (1 - β_ij²).

The aggregate clause `‖Vᵀ V̂_svdstack‖_F² → Σ_j S_j/(S_j+1)` is a corollary of the E2 form of
`thm:gen_rank_weight_svdstak`, which delivers the limit as `tr(B_Rᵀ A_{β,R}⁻¹ B_R)`. This file
proves the deterministic identity that connects the two:

    tr(B_Rᵀ A_{β,R}⁻¹ B_R) = Σ_j S_j / (S_j + 1)   in the exactly aligned model.

"Exactly aligned" means `r_i = r` for every table and `R_i = 1`, so the `j`-th spike of every
table is the shared direction `v_j`. In the block objects of `RankR/General.lean` this is
`rk := alignedRk M r` and `R := alignedR M r`.

## Route

1. `inv_add_mul_transpose_mul`: the push-through identity
   `(D + B Bᵀ)⁻¹ B = D⁻¹ B (1 + Bᵀ D⁻¹ B)⁻¹`, valid for every positive definite `D` and every
   `B`. It follows from `(D + B Bᵀ) D⁻¹ B (1 + G)⁻¹ = (B + B G)(1 + G)⁻¹ = B` with
   `G = Bᵀ D⁻¹ B`.
2. `trace_conj_inv_add_eq`: the trace form `tr(Bᵀ (D + B Bᵀ)⁻¹ B) = tr(G (1 + G)⁻¹)`.
3. `transpose_mul_inv_DBlock_mul_aligned`: in the aligned model `G = diag(S_1, …, S_r)`.
4. `trace_conj_inv_ABlock_aligned`: the two combine to `Σ_j S_j/(S_j+1)`.
5. `trace_conj_inv_AbetaR_aligned`: the same statement on the `r_i = 1` objects `BR`, `AbetaR`
   of `RankR/Defs.lean` through the `rfl` bridges of `RankR/Flatten.lean`, which is the shape
   the E2 form of `thm:gen_rank_weight_svdstak` reads.

Everything here is deterministic linear algebra. No probability and no random matrix theory
limit law enters.

STATUS 2026-09-01: proved, 0 `sorry`.
-/

open scoped Matrix

namespace StackedSVD

/-! ### 1. The push-through identity

`D ≻ 0` and `B` arbitrary. Nothing in this section knows about the block model. -/

section PushThrough

variable {p q : ℕ} {D : Matrix (Fin p) (Fin p) ℝ}

/-- `1 + Bᵀ D⁻¹ B ≻ 0` when `D ≻ 0`, so the middle factor of the push-through identity is
invertible. -/
theorem posDef_one_add_conj_inv (hD : D.PosDef) (B : Matrix (Fin p) (Fin q) ℝ) :
    ((1 : Matrix (Fin q) (Fin q) ℝ) + Bᵀ * D⁻¹ * B).PosDef := by
  have hpsd : (Bᵀ * D⁻¹ * B).PosSemidef := by
    have h := hD.inv.posSemidef.conjTranspose_mul_mul_same B
    rwa [Matrix.conjTranspose_eq_transpose_of_trivial] at h
  rw [add_comm]
  exact posSemidef_add_posDef hpsd Matrix.PosDef.one

/-- `D + B Bᵀ ≻ 0` when `D ≻ 0`. -/
theorem posDef_add_mul_transpose (hD : D.PosDef) (B : Matrix (Fin p) (Fin q) ℝ) :
    (D + B * Bᵀ).PosDef := by
  rw [add_comm]
  exact posSemidef_add_posDef
    (by simpa using Matrix.posSemidef_self_mul_conjTranspose B) hD

/-- **Push-through identity.** `(D + B Bᵀ)⁻¹ B = D⁻¹ B (1 + Bᵀ D⁻¹ B)⁻¹` for every positive
definite `D`. This is the Woodbury step of `paper_edits.md` E4 item 1. -/
theorem inv_add_mul_transpose_mul (hD : D.PosDef) (B : Matrix (Fin p) (Fin q) ℝ) :
    (D + B * Bᵀ)⁻¹ * B = D⁻¹ * B * (1 + Bᵀ * D⁻¹ * B)⁻¹ := by
  have hA := posDef_add_mul_transpose hD B
  have hG := posDef_one_add_conj_inv hD B
  have hDdet : IsUnit D.det := (Matrix.isUnit_iff_isUnit_det D).mp hD.isUnit
  have hAdet : IsUnit (D + B * Bᵀ).det := (Matrix.isUnit_iff_isUnit_det _).mp hA.isUnit
  have hGdet : IsUnit ((1 : Matrix (Fin q) (Fin q) ℝ) + Bᵀ * D⁻¹ * B).det :=
    (Matrix.isUnit_iff_isUnit_det _).mp hG.isUnit
  -- `(D + B Bᵀ) D⁻¹ B = B (1 + Bᵀ D⁻¹ B)`
  have h1 : (D + B * Bᵀ) * (D⁻¹ * B) = B * (1 + Bᵀ * D⁻¹ * B) := by
    rw [Matrix.add_mul, ← Matrix.mul_assoc D, Matrix.mul_nonsing_inv D hDdet, Matrix.one_mul,
      Matrix.mul_add, Matrix.mul_one]
    simp [Matrix.mul_assoc]
  -- so `(D + B Bᵀ) (D⁻¹ B (1 + G)⁻¹) = B`
  have key : (D + B * Bᵀ) * (D⁻¹ * B * (1 + Bᵀ * D⁻¹ * B)⁻¹) = B := by
    calc (D + B * Bᵀ) * (D⁻¹ * B * (1 + Bᵀ * D⁻¹ * B)⁻¹)
        = ((D + B * Bᵀ) * (D⁻¹ * B)) * (1 + Bᵀ * D⁻¹ * B)⁻¹ := by rw [← Matrix.mul_assoc]
      _ = B * ((1 + Bᵀ * D⁻¹ * B) * (1 + Bᵀ * D⁻¹ * B)⁻¹) := by rw [h1, Matrix.mul_assoc]
      _ = B := by rw [Matrix.mul_nonsing_inv _ hGdet, Matrix.mul_one]
  have h2 : (D + B * Bᵀ)⁻¹ * ((D + B * Bᵀ) * (D⁻¹ * B * (1 + Bᵀ * D⁻¹ * B)⁻¹))
      = (D + B * Bᵀ)⁻¹ * B := by rw [key]
  rw [← Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hAdet, Matrix.one_mul] at h2
  exact h2.symm

/-- The trace form of the push-through identity:
`tr(Bᵀ (D + B Bᵀ)⁻¹ B) = tr(G (1 + G)⁻¹)` with `G = Bᵀ D⁻¹ B`. -/
theorem trace_conj_inv_add_eq (hD : D.PosDef) (B : Matrix (Fin p) (Fin q) ℝ) :
    Matrix.trace (Bᵀ * (D + B * Bᵀ)⁻¹ * B)
      = Matrix.trace ((Bᵀ * D⁻¹ * B) * (1 + Bᵀ * D⁻¹ * B)⁻¹) := by
  rw [Matrix.mul_assoc, inv_add_mul_transpose_mul hD B]
  congr 1
  simp [Matrix.mul_assoc]

end PushThrough

/-! ### 2. Diagonal inverses and the flat sum

Two small tools. `inv_diagonal_of_ne_zero` is the field form of `Matrix.inv_diagonal`, proved
through `Matrix.inv_eq_right_inv` so that no `Ring.inverse` appears. `sum_flat` turns a sum
over the flat index `Fin (r̃)` into the double sum over `(i, j)` that the paper writes. -/

/-- The inverse of a diagonal matrix with nonzero entries. -/
theorem inv_diagonal_of_ne_zero {p : ℕ} (v : Fin p → ℝ) (hv : ∀ i, v i ≠ 0) :
    (Matrix.diagonal v)⁻¹ = Matrix.diagonal fun i => (v i)⁻¹ := by
  refine Matrix.inv_eq_right_inv ?_
  rw [Matrix.diagonal_mul_diagonal,
    show (fun i => v i * (v i)⁻¹) = fun _ : Fin p => (1 : ℝ) from
      funext fun i => mul_inv_cancel₀ (hv i)]
  exact Matrix.diagonal_one

variable {M r : ℕ} {rk : Fin M → ℕ}

/-- A sum over the flat index `p ∈ [r̃]` is the double sum over the block index `(i, j)`. -/
theorem sum_flat {γ : Type*} [AddCommMonoid γ] (f : Fin (rtot rk) → γ) :
    ∑ p, f p = ∑ i, ∑ j, f (flat i j) := by
  rw [← Fintype.sum_equiv finSigmaFinEquiv (fun s => f (finSigmaFinEquiv s)) f fun _ => rfl,
    Fintype.sum_sigma]
  rfl

/-- `β` read at a flat index built from a block index. -/
theorem betaFlat_flat (β : (i : Fin M) → Fin (rk i) → ℝ) (i : Fin M) (j : Fin (rk i)) :
    betaFlat β (flat i j) = β i j :=
  congrArg (fun s : (i : Fin M) × Fin (rk i) => β s.1 s.2) (blk_flat i j)

/-- `D⁻¹ = diag((1 - β_ij²)⁻¹)` whenever every `β_ij² < 1`. -/
theorem inv_DBlock {β : (i : Fin M) → Fin (rk i) → ℝ} (h1 : ∀ i j, β i j ^ 2 < 1) :
    (DBlock (rk := rk) β)⁻¹ = Matrix.diagonal fun p => (1 - betaFlat β p ^ 2)⁻¹ := by
  refine inv_diagonal_of_ne_zero _ fun p => ?_
  have h := h1 (blk p).1 (blk p).2
  rw [betaFlat]
  intro hz
  rw [sub_eq_zero] at hz
  exact absurd hz.symm (ne_of_lt h)

/-! ### 3. The exactly aligned family

`r_i = r` for every table and `R_i = 1`: the `j`-th spike of table `i` is the shared direction
`v_j`. This is the model of `thm:rank_r_svdstack` (`main_paper.tex:2391`). -/

/-- Every table carries `r` spikes. -/
abbrev alignedRk (M r : ℕ) : Fin M → ℕ := fun _ => r

/-- Every alignment matrix is the identity, so spike `j` of table `i` is the shared `v_j`. -/
def alignedR (M r : ℕ) : (i : Fin M) → Matrix (Fin r) (Fin (alignedRk M r i)) ℝ := fun _ => 1

/-- `S_j = Σ_i β_ij² / (1 - β_ij²)`, the paper's aggregate signal strength of component `j`
(`main_paper.tex:2395`). At `r = 1` this is `Sval` of `SVDStack/Defs.lean`. -/
noncomputable def Sagg (β : Fin M → Fin r → ℝ) (j : Fin r) : ℝ :=
  ∑ i, β i j ^ 2 / (1 - β i j ^ 2)

/-- `S_j ≥ 0` when every `β_ij² < 1`. -/
theorem Sagg_nonneg {β : Fin M → Fin r → ℝ} (h1 : ∀ i j, β i j ^ 2 < 1) (j : Fin r) :
    0 ≤ Sagg β j :=
  Finset.sum_nonneg fun i _ => div_nonneg (sq_nonneg _) (by linarith [h1 i j])

/-- Entry of `B_R` in the aligned model: row `(i, j)` has the single nonzero entry `β_ij` in
column `j`. -/
theorem BBlock_aligned_apply (β : Fin M → Fin r → ℝ) (p : Fin (rtot (alignedRk M r)))
    (k : Fin r) :
    BBlock (rk := alignedRk M r) β (alignedR M r) p k =
      if k = (blk p).2 then betaFlat (rk := alignedRk M r) β p else 0 := by
  simp only [BBlock, Matrix.of_apply, alignedR, Matrix.one_apply]
  split_ifs <;> ring

/-- Entry of `B_R` at a flat index built from a block index, in the aligned model. Row `(i, j)`
carries the single nonzero entry `β_ij` in column `j`. -/
theorem BBlock_aligned_flat (β : Fin M → Fin r → ℝ) (i : Fin M) (j k : Fin r) :
    BBlock (rk := alignedRk M r) β (alignedR M r) (flat i j) k = if k = j then β i j else 0 := by
  rw [BBlock_aligned_apply]
  exact congrArg (fun s : (i : Fin M) × Fin (alignedRk M r i) =>
    (if k = (s.2 : Fin r) then β s.1 (s.2 : Fin r) else 0 : ℝ))
    (blk_flat (rk := alignedRk M r) i j)

/-- **In the aligned model `Bᵀ D⁻¹ B = diag(S_1, …, S_r)`.** -/
theorem transpose_mul_inv_DBlock_mul_aligned (β : Fin M → Fin r → ℝ)
    (h1 : ∀ i j, β i j ^ 2 < 1) :
    (BBlock (rk := alignedRk M r) β (alignedR M r))ᵀ * (DBlock (rk := alignedRk M r) β)⁻¹ *
        BBlock (rk := alignedRk M r) β (alignedR M r) = Matrix.diagonal (Sagg β) := by
  have h1' : ∀ (i : Fin M) (j : Fin (alignedRk M r i)), β i j ^ 2 < 1 := h1
  ext k l
  rw [inv_DBlock h1', Matrix.diagonal_apply, Matrix.mul_apply, sum_flat (rk := alignedRk M r)]
  simp only [Matrix.mul_diagonal, Matrix.transpose_apply, BBlock_aligned_flat, betaFlat_flat]
  by_cases hkl : k = l
  · subst hkl
    rw [if_pos rfl, Sagg]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Finset.sum_congr rfl fun j (_ : j ∈ (Finset.univ : Finset (Fin r))) =>
        (show (if k = j then β i j else 0) * (1 - β i j ^ 2)⁻¹ * (if k = j then β i j else 0)
            = if k = j then β i j ^ 2 * (1 - β i j ^ 2)⁻¹ else 0 by split_ifs <;> ring),
      Finset.sum_ite_eq Finset.univ k fun j => β i j ^ 2 * (1 - β i j ^ 2)⁻¹,
      if_pos (Finset.mem_univ k), div_eq_mul_inv]
  · rw [if_neg hkl]
    refine Finset.sum_eq_zero fun i _ => Finset.sum_eq_zero fun j _ => ?_
    by_cases hkj : k = j
    · have hlj : l ≠ j := fun h => hkl (hkj.trans h.symm)
      rw [if_neg hlj, mul_zero]
    · rw [if_neg hkj, zero_mul, zero_mul]

/-- **The block-diagonal Woodbury identity** (`paper_edits.md` E4 item 1). In the exactly
aligned model of `thm:rank_r_svdstack` (`main_paper.tex:2391`),
`tr(B_Rᵀ A_{β,R}⁻¹ B_R) = Σ_j S_j / (S_j + 1)`. -/
theorem trace_conj_inv_ABlock_aligned (β : Fin M → Fin r → ℝ) (h0 : ∀ i j, 0 ≤ β i j)
    (h1 : ∀ i j, β i j < 1) :
    Matrix.trace ((BBlock (rk := alignedRk M r) β (alignedR M r))ᵀ *
        (ABlock (rk := alignedRk M r) β (alignedR M r))⁻¹ *
        BBlock (rk := alignedRk M r) β (alignedR M r))
      = ∑ j, Sagg β j / (Sagg β j + 1) := by
  have h0' : ∀ (i : Fin M) (j : Fin (alignedRk M r i)), 0 ≤ β i j := h0
  have h1' : ∀ (i : Fin M) (j : Fin (alignedRk M r i)), β i j < 1 := h1
  have hsq : ∀ i j, β i j ^ 2 < 1 := fun i j => by nlinarith [h0 i j, h1 i j]
  have hD : (DBlock (rk := alignedRk M r) β).PosDef := DBlock_posDef h0' h1'
  have hA : ABlock (rk := alignedRk M r) β (alignedR M r)
      = DBlock (rk := alignedRk M r) β +
          BBlock (rk := alignedRk M r) β (alignedR M r) *
            (BBlock (rk := alignedRk M r) β (alignedR M r))ᵀ := by
    rw [ABlock, add_comm]
  rw [hA, trace_conj_inv_add_eq hD, transpose_mul_inv_DBlock_mul_aligned β hsq]
  have hone : (1 : Matrix (Fin r) (Fin r) ℝ) + Matrix.diagonal (Sagg β)
      = Matrix.diagonal fun j => Sagg β j + 1 := by
    rw [← Matrix.diagonal_one, Matrix.diagonal_add]
    exact congrArg Matrix.diagonal (funext fun j => add_comm _ _)
  rw [hone, inv_diagonal_of_ne_zero _ fun j => by
      have := Sagg_nonneg hsq j; positivity,
    Matrix.diagonal_mul_diagonal, Matrix.trace_diagonal]
  exact Finset.sum_congr rfl fun j _ => (div_eq_mul_inv _ _).symm

/-- The same identity on the `r_i = 1` objects `BR`, `AbetaR` of `RankR/Defs.lean`, through the
`rfl` bridges `BBlock_eq_BR` and `ABlock_eq_AbetaR` of `RankR/Flatten.lean`. The E2 form of
`thm:gen_rank_weight_svdstak` is stated in these objects and reads this directly. -/
theorem trace_conj_inv_AbetaR_aligned (β : Fin M → Fin r → ℝ) (h0 : ∀ i j, 0 ≤ β i j)
    (h1 : ∀ i j, β i j < 1) :
    Matrix.trace ((BR (betaFlat (rk := alignedRk M r) β) (Rcol (alignedR M r)))ᵀ *
        (AbetaR (betaFlat (rk := alignedRk M r) β) (Rcol (alignedR M r)))⁻¹ *
        BR (betaFlat (rk := alignedRk M r) β) (Rcol (alignedR M r)))
      = ∑ j, Sagg β j / (Sagg β j + 1) :=
  trace_conj_inv_ABlock_aligned β h0 h1

/-! ### 4. `S_j` in the paper's data `(θ, c)`

`β_ij² = betaSq θ_ij c_i` (`Defs.lean:41`), so `S_j = 0` exactly when every table is
subcritical at component `j`. This is the case `paper_edits.md` E4 item 2 excludes from the
component-labelled clause. -/

/-- `betaSq θ c = 0` exactly below the threshold, when `c > 0`. -/
theorem betaSq_eq_zero_iff {θ c : ℝ} (hc : 0 < c) : betaSq θ c = 0 ↔ ¬ (θ ^ 4 > c) := by
  constructor
  · intro h hgt
    rw [betaSq, if_pos hgt] at h
    have hd : (0 : ℝ) < θ ^ 4 + θ ^ 2 := by nlinarith [sq_nonneg θ, sq_nonneg (θ ^ 2)]
    rw [div_eq_zero_iff] at h
    rcases h with h | h
    · linarith
    · linarith
  · intro h
    rw [betaSq, if_neg h]

/-- `S_j = 0` exactly when every table is subcritical at component `j`, that is
`θ_ij⁴ ≤ c_i` for every `i`. -/
theorem Sagg_eq_zero_iff (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (j : Fin r) :
    Sagg (fun i j => Real.sqrt (betaSq (θ i j) (c i))) j = 0 ↔ ∀ i, ¬ (θ i j ^ 4 > c i) := by
  have hsq : ∀ i : Fin M, Real.sqrt (betaSq (θ i j) (c i)) ^ 2 = betaSq (θ i j) (c i) :=
    fun i => Real.sq_sqrt (Scalars.betaSq_nonneg _ _)
  have hpos : ∀ i : Fin M, 0 < 1 - betaSq (θ i j) (c i) :=
    fun i => Scalars.one_sub_betaSq_pos (hc i)
  have hterm : ∀ i : Fin M,
      Real.sqrt (betaSq (θ i j) (c i)) ^ 2 / (1 - Real.sqrt (betaSq (θ i j) (c i)) ^ 2)
        = betaSq (θ i j) (c i) / (1 - betaSq (θ i j) (c i)) := by
    intro i; rw [hsq i]
  rw [Sagg, Finset.sum_congr rfl fun i (_ : i ∈ Finset.univ) => hterm i,
    Finset.sum_eq_zero_iff_of_nonneg fun i _ =>
      div_nonneg (Scalars.betaSq_nonneg _ _) (hpos i).le]
  constructor
  · intro h i
    have hi := h i (Finset.mem_univ i)
    rw [div_eq_zero_iff] at hi
    rcases hi with hi | hi
    · exact (betaSq_eq_zero_iff (hc i)).mp hi
    · linarith [hpos i]
  · intro h i _
    rw [(betaSq_eq_zero_iff (hc i)).mpr (h i), zero_div]

end StackedSVD
