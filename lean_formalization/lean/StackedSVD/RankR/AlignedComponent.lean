/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.Aligned
import StackedSVD.RankR.WeightedUpperG
import StackedSVD.RankR.Flatten
import StackedSVD.LinAlg.SpecIdx
import StackedSVD.LinAlg.SpecProjPerturb
import StackedSVD.LinAlg.Eigen

/-!
# The aligned limit spectrum at one component index

Step 2 of `notes/archive/rankr_D4_plan.md` (items 2.1, 2.2, 2.4 to 2.7 and 6.3), the deterministic
half of the component clause of `thm:rank_r_svdstack` (`main_paper.tex:2405`). Everything here
is linear algebra: no probability and no random matrix theory limit law enters.

Write `W⋆ = optWG β`, `A⋆ = W⋆ A_{β,R} W⋆ᵀ`, `C = W⋆ B_R` in the exactly aligned model
(`r_i = r`, `R_i = I`), and `S_j = Sagg β j`. The chain is

1. `ABlockW_optWG_aligned_eq`: `A⋆ = C Cᵀ + I` (`main_paper.tex:2103`).
2. `transpose_BBlockW_optWG_aligned`: `Cᵀ C = diag(S_0, …, S_{r-1})`, so the columns of `C`
   are orthogonal with `‖c_j‖² = S_j`.
3. `eigenvalues₀_ABlockW_optWG_aligned`: the sorted spectrum of `A⋆` is
   `1 + S_0, …, 1 + S_{r-1}, 1, …, 1`.
4. `topGap_ABlockW_optWG_aligned` and its successor, and `simpleSpec_ABlockW_optWG_aligned`:
   under `StrictAnti (Sagg β)` and `0 < S_j` the index `j` is separated on both sides and
   every index at or below `j` is simple. Each of these, and items 5 and 6 below, has an
   `_of_sep` twin that takes only `SaggSep β`, the separation that `Sagg` really needs:
   `Antitone (Sagg β)` plus a strict drop from every component with `S > 0`. The `StrictAnti`
   forms are one-line wrappers (`Sagg_sep_of_strictAnti`), and section 4 proves `SaggSep`
   inside the model with no hypothesis (b) at all (`UnalignedModelR.saggSep_of_model`).
5. `normSq_specProjIdx_ABlockW_aligned`: `‖P_j y‖² = ⟪c_j, y⟫² / S_j`, where `P_j` is the
   projector at the sorted eigenvalue index `j` (`LinAlg/SpecIdx.lean`).
6. `compLimit_aligned`: the value the component clause converges to,
   `‖P_j c_j‖² / (1 + S_j) = S_j / (S_j + 1)`.
7. Item 6.3, section 4: `Sagg_antitone` and `Sagg_lt_of_pos`, with the model-facing pair
   `UnalignedModelR.Sagg_antitone_of_model` and `Sagg_lt_of_pos_of_model`. They read the
   spike order `SpikedModelR.hθanti` and the monotonicity of `θ ↦ (θ⁴-c)/(θ²+c)`, which is
   finding F1 of the plan: the paper's hypothesis (b) `S_j ≠ S_k` follows from the sorted
   form of hypothesis (a) that the Lean model carries (`hθanti`, the same strict order in
   every table) at every index where the component clause has content. At unordered `θ`
   (a) alone does not give (b): plan cell D has `θ` unordered and `S_0 = S_1`.

Modeling choice M5 of the plan: `A⋆` is **not** block diagonal in the flat index, which groups
by table first. The orthogonal-column form `A⋆ = I + ∑_j c_j c_jᵀ` needs no permutation and is
one rewrite from `transpose_mul_inv_DBlock_mul_aligned` (`RankR/Aligned.lean`).

Sections 1 to 3 take the deterministic hypotheses `h0 : ∀ i j, 0 ≤ β i j` and
`h1 : ∀ i j, β i j < 1` of `RankR/Aligned.lean`, not the model-facing `hc` and `hβdef`: no
statement there mentions a model. `UnalignedModelR.thm_rank_r_svdstack_aggregate`
(`RankR/WeightedUpperG.lean:638`) derives `h0` and `h1` from `hβdef` and `hc` in two lines.
Only the two corollaries at the end of section 4 read a model.

STATUS 2026-09-02: proved, 0 `sorry`.
-/

open Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

variable {M r : ℕ}

/-! ### 1. The two matrix identities at `W⋆` -/

section Identities

variable {rk : Fin M → ℕ}

/-- Column `k` of `W B_R`, as a vector of `EuclideanSpace`. Mirror: `Rcol`
(`RankR/Flatten.lean`) and the planned `VtVWGcol` of the probabilistic layer. -/
noncomputable def BBlockWcol (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ) (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ)
    (k : Fin r) : EuclideanSpace ℝ (Fin (rtot rk)) :=
  WithLp.toLp 2 fun p => BBlockW W β R p k

theorem ofLp_BBlockWcol (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ) (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ)
    (k : Fin r) :
    WithLp.ofLp (BBlockWcol W β R k) = BBlockW W β R *ᵥ Pi.single k 1 := by
  rw [Matrix.mulVec_single_one]
  rfl

/-- `W⋆ᵀ W⋆ = D⁻¹`: the optimal weight is `D^{-1/2}` and is diagonal, hence symmetric. -/
theorem transpose_optWG_mul_optWG {β : (i : Fin M) → Fin (rk i) → ℝ}
    (h1 : ∀ i j, β i j ^ 2 < 1) :
    (optWG (rk := rk) β)ᵀ * optWG (rk := rk) β = (DBlock (rk := rk) β)⁻¹ := by
  have h1' : ∀ p, betaFlat (rk := rk) β p ^ 2 < 1 := fun p => h1 (blk p).1 (blk p).2
  rw [optWG, Matrix.diagonal_transpose, Matrix.diagonal_mul_diagonal, inv_DBlock h1]
  refine congrArg Matrix.diagonal (funext fun p => ?_)
  have hne : (1 : ℝ) - betaFlat (rk := rk) β p ^ 2 ≠ 0 := by
    have := h1' p; intro hz; rw [sub_eq_zero] at hz; linarith [hz.symm]
  have hkey := optW_mul_one_sub_mul (β := betaFlat (rk := rk) β) (i := p) (h1' p)
  field_simp
  nlinarith [hkey]

/-- **Item 2.1**: `A⋆ = C Cᵀ + I` in the exactly aligned model (`main_paper.tex:2103`), read
through the flatten bridges of `RankR/Flatten.lean`. -/
theorem ABlockW_optWG_aligned_eq (β : Fin M → Fin r → ℝ) (h1 : ∀ i j, β i j ^ 2 < 1) :
    ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r)
      = BBlockW (optWG (rk := alignedRk M r) β) β (alignedR M r) *
          (BBlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))ᵀ + 1 := by
  rw [ABlockW_eq_AbetaRW, BBlockW_eq_BRW, optWG_eq_optWR]
  exact abetaRW_optWR_eq (Rcol (alignedR M r)) fun p => h1 (blk p).1 (blk p).2

/-- **Item 2.2**: `Cᵀ C = diag(S_0, …, S_{r-1})` in the exactly aligned model. Route:
`C = W⋆ B_R`, `W⋆ᵀ W⋆ = D⁻¹`, then `transpose_mul_inv_DBlock_mul_aligned`
(`RankR/Aligned.lean`). -/
theorem transpose_BBlockW_optWG_aligned (β : Fin M → Fin r → ℝ) (h1 : ∀ i j, β i j ^ 2 < 1) :
    (BBlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))ᵀ *
        BBlockW (optWG (rk := alignedRk M r) β) β (alignedR M r)
      = Matrix.diagonal (Sagg β) := by
  have hkey : (BBlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))ᵀ *
      BBlockW (optWG (rk := alignedRk M r) β) β (alignedR M r)
      = (BBlock (rk := alignedRk M r) β (alignedR M r))ᵀ *
          ((optWG (rk := alignedRk M r) β)ᵀ * optWG (rk := alignedRk M r) β) *
          BBlock (rk := alignedRk M r) β (alignedR M r) := by
    simp only [BBlockW, Matrix.transpose_mul, Matrix.mul_assoc]
  rw [hkey, transpose_optWG_mul_optWG h1, transpose_mul_inv_DBlock_mul_aligned β h1]

end Identities

/-! ### 2. The sorted spectrum of `A⋆` -/

section Spectrum

/-- **Item 2.4**: the sorted spectrum of `A⋆` is `1 + S` padded by ones. Route: item 2.1,
`eigenvalues₀_add_one` and `eigenvalues₀_mul_transpose_of_transpose_mul_diagonal`
(`LinAlg/Eigen.lean`) on top of item 2.2. -/
theorem eigenvalues₀_ABlockW_optWG_aligned (β : Fin M → Fin r → ℝ) (h0 : ∀ i j, 0 ≤ β i j)
    (h1 : ∀ i j, β i j < 1) (hanti : Antitone (Sagg β)) (hM : 0 < M)
    (k : Fin (Fintype.card (Fin (rtot (alignedRk M r))))) :
    (isHermitian_ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r)).eigenvalues₀ k
      = 1 + (if h : (k : ℕ) < r then Sagg β ⟨(k : ℕ), h⟩ else 0) := by
  have hsq : ∀ i j, β i j ^ 2 < 1 := fun i j => by nlinarith [h0 i j, h1 i j]
  set C := BBlockW (optWG (rk := alignedRk M r) β) β (alignedR M r) with hC
  have hP : (C * Cᵀ).IsHermitian := isHermitian_mul_transpose_self C
  have hP1 : (C * Cᵀ + 1).IsHermitian := hP.add Matrix.isHermitian_one
  have heq : ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r) = C * Cᵀ + 1 :=
    ABlockW_optWG_aligned_eq β hsq
  have hcongr := eigenvalues₀_congr_mat
    (isHermitian_ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r)) hP1 heq
  rw [congrFun hcongr k, eigenvalues₀_add_one hP hP1 k,
    eigenvalues₀_mul_transpose_of_transpose_mul_diagonal C
      (transpose_BBlockW_optWG_aligned β hsq) hanti (Sagg_nonneg hsq)
      (le_rtot_alignedRk hM) hP k,
    add_comm]

/-- The sorted index `j` of a component lies inside the flat index range. -/
theorem lt_card_rtot_alignedRk (hM : 0 < M) {j : ℕ} (hj : j < r) :
    j < Fintype.card (Fin (rtot (alignedRk M r))) := by
  simpa using lt_of_lt_of_le hj (le_rtot_alignedRk (M := M) (r := r) hM)

/-- The separation `Sagg` really needs: `S` decreases strictly from every component that
carries `S > 0`. It is what `StrictAnti (Sagg β)` gives (`Sagg_sep_of_strictAnti`) and what
`Sagg_lt_of_pos` of section 4 proves inside the model with no `hS` at all. Every lemma of
item 2.5 goes through this form, so the `hS`-free facade of the plan's step 6.3 needs no edit
to this file. -/
def SaggSep (β : Fin M → Fin r → ℝ) : Prop :=
  Antitone (Sagg β) ∧ ∀ a b : Fin r, a < b → 0 < Sagg β a → Sagg β b < Sagg β a

theorem Sagg_sep_of_strictAnti {β : Fin M → Fin r → ℝ} (hS : StrictAnti (Sagg β)) :
    SaggSep β := ⟨hS.antitone, fun _ _ hab _ => hS hab⟩

/-- The eigenvalue at index `b` is strictly below the one at index `a` whenever `a < b` and
`a` is at or before a component with `S_a > 0`. The single monotonicity fact behind the two
gaps and the simplicity of item 2.5. -/
theorem padSagg_lt {β : Fin M → Fin r → ℝ} (hS : SaggSep β) {j : Fin r}
    (hpos : 0 < Sagg β j) {a b : ℕ} (hab : a < b) (haj : a ≤ (j : ℕ)) :
    1 + (if h : b < r then Sagg β ⟨b, h⟩ else 0)
      < 1 + (if h : a < r then Sagg β ⟨a, h⟩ else 0) := by
  have har : a < r := lt_of_le_of_lt haj j.isLt
  have hSa : 0 < Sagg β ⟨a, har⟩ :=
    lt_of_lt_of_le hpos (hS.1 (show (⟨a, har⟩ : Fin r) ≤ j from haj))
  rw [dif_pos har]
  by_cases hb : b < r
  · rw [dif_pos hb]
    have hlt : Sagg β ⟨b, hb⟩ < Sagg β ⟨a, har⟩ :=
      hS.2 ⟨a, har⟩ ⟨b, hb⟩ hab hSa
    linarith
  · rw [dif_neg hb]
    linarith

/-- **Item 2.5, the gap below index `j`**: `TopGap A⋆ j`. -/
theorem topGap_ABlockW_optWG_aligned_of_sep (β : Fin M → Fin r → ℝ) (h0 : ∀ i j, 0 ≤ β i j)
    (h1 : ∀ i j, β i j < 1) (hsep : SaggSep β) (hM : 0 < M) (j : Fin r)
    (hpos : 0 < Sagg β j) :
    TopGap (ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))
      (isHermitian_ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r)) (j : ℕ) := by
  intro k l hk hl
  rw [eigenvalues₀_ABlockW_optWG_aligned β h0 h1 hsep.1 hM,
    eigenvalues₀_ABlockW_optWG_aligned β h0 h1 hsep.1 hM]
  exact padSagg_lt hsep hpos (lt_of_lt_of_le hk hl) (le_of_lt hk)

/-- **Item 2.5, the gap below index `j`** at the paper's hypothesis (b). -/
theorem topGap_ABlockW_optWG_aligned (β : Fin M → Fin r → ℝ) (h0 : ∀ i j, 0 ≤ β i j)
    (h1 : ∀ i j, β i j < 1) (hS : StrictAnti (Sagg β)) (hM : 0 < M) (j : Fin r)
    (hpos : 0 < Sagg β j) :
    TopGap (ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))
      (isHermitian_ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r)) (j : ℕ) :=
  topGap_ABlockW_optWG_aligned_of_sep β h0 h1 (Sagg_sep_of_strictAnti hS) hM j hpos

/-- **Item 2.5, the gap above index `j`**: `TopGap A⋆ (j + 1)`. -/
theorem topGap_succ_ABlockW_optWG_aligned_of_sep (β : Fin M → Fin r → ℝ)
    (h0 : ∀ i j, 0 ≤ β i j) (h1 : ∀ i j, β i j < 1) (hsep : SaggSep β) (hM : 0 < M) (j : Fin r)
    (hpos : 0 < Sagg β j) :
    TopGap (ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))
      (isHermitian_ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r)) ((j : ℕ) + 1) := by
  intro k l hk hl
  rw [eigenvalues₀_ABlockW_optWG_aligned β h0 h1 hsep.1 hM,
    eigenvalues₀_ABlockW_optWG_aligned β h0 h1 hsep.1 hM]
  exact padSagg_lt hsep hpos (lt_of_lt_of_le hk hl) (by omega)

/-- **Item 2.5, the gap above index `j`** at the paper's hypothesis (b). -/
theorem topGap_succ_ABlockW_optWG_aligned (β : Fin M → Fin r → ℝ) (h0 : ∀ i j, 0 ≤ β i j)
    (h1 : ∀ i j, β i j < 1) (hS : StrictAnti (Sagg β)) (hM : 0 < M) (j : Fin r)
    (hpos : 0 < Sagg β j) :
    TopGap (ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))
      (isHermitian_ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r)) ((j : ℕ) + 1) :=
  topGap_succ_ABlockW_optWG_aligned_of_sep β h0 h1 (Sagg_sep_of_strictAnti hS) hM j hpos

/-- **Item 2.5, simplicity**: every sorted index at or below `j` carries a simple eigenvalue
of `A⋆`. This is what makes `specProjIdx A⋆ j` a rank-one projector. -/
theorem simpleSpec_ABlockW_optWG_aligned_of_sep (β : Fin M → Fin r → ℝ)
    (h0 : ∀ i j, 0 ≤ β i j) (h1 : ∀ i j, β i j < 1) (hsep : SaggSep β) (hM : 0 < M) (j : Fin r)
    (hpos : 0 < Sagg β j) :
    SimpleSpec (ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))
      (isHermitian_ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))
      ((j : ℕ) + 1) := by
  intro k l hk hkl
  rw [eigenvalues₀_ABlockW_optWG_aligned β h0 h1 hsep.1 hM,
    eigenvalues₀_ABlockW_optWG_aligned β h0 h1 hsep.1 hM]
  have hne : (k : ℕ) ≠ (l : ℕ) := fun h => hkl (Fin.val_injective h)
  rcases lt_or_gt_of_ne hne with h | h
  · exact ne_of_lt (padSagg_lt hsep hpos h (by omega))
  · exact ne_of_gt (padSagg_lt hsep hpos h (by omega))

/-- **Item 2.5, simplicity** at the paper's hypothesis (b). -/
theorem simpleSpec_ABlockW_optWG_aligned (β : Fin M → Fin r → ℝ) (h0 : ∀ i j, 0 ≤ β i j)
    (h1 : ∀ i j, β i j < 1) (hS : StrictAnti (Sagg β)) (hM : 0 < M) (j : Fin r)
    (hpos : 0 < Sagg β j) :
    SimpleSpec (ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))
      (isHermitian_ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))
      ((j : ℕ) + 1) :=
  simpleSpec_ABlockW_optWG_aligned_of_sep β h0 h1 (Sagg_sep_of_strictAnti hS) hM j hpos

/-- **Item 2.5, the eigenvalue set at index `j`** is the singleton `{1 + S_j}`. -/
theorem eigSetIdx_ABlockW_optWG_aligned (β : Fin M → Fin r → ℝ) (h0 : ∀ i j, 0 ≤ β i j)
    (h1 : ∀ i j, β i j < 1) (hanti : Antitone (Sagg β)) (hM : 0 < M) (j : Fin r) :
    eigSetIdx (ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))
        (isHermitian_ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r)) (j : ℕ)
      = {1 + Sagg β j} := by
  have hjc : (j : ℕ) < Fintype.card (Fin (rtot (alignedRk M r))) :=
    lt_card_rtot_alignedRk hM j.isLt
  rw [eigSetIdx_eq_singleton _ _ hjc,
    eigenvalues₀_ABlockW_optWG_aligned β h0 h1 hanti hM ⟨(j : ℕ), hjc⟩, dif_pos j.isLt]

end Spectrum

/-! ### 3. The index projector at a supercritical component -/

section Projector

/-- `‖c_j‖² = S_j`: the columns of `C = W⋆ B_R` have squared norm `S_j` (item 2.2 on the
diagonal). -/
theorem inner_self_BBlockWcol_aligned (β : Fin M → Fin r → ℝ) (h1 : ∀ i j, β i j ^ 2 < 1)
    (j : Fin r) :
    ⟪BBlockWcol (optWG (rk := alignedRk M r) β) β (alignedR M r) j,
        BBlockWcol (optWG (rk := alignedRk M r) β) β (alignedR M r) j⟫_ℝ = Sagg β j := by
  have hCC := transpose_BBlockW_optWG_aligned β h1
  have hentry : ((BBlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))ᵀ *
      BBlockW (optWG (rk := alignedRk M r) β) β (alignedR M r)) j j = Sagg β j := by
    rw [hCC, Matrix.diagonal_apply_eq]
  rw [real_inner_eq_dotProduct, ← hentry, Matrix.mul_apply]
  exact Finset.sum_congr rfl fun p _ => rfl

theorem norm_sq_BBlockWcol_aligned (β : Fin M → Fin r → ℝ) (h1 : ∀ i j, β i j ^ 2 < 1)
    (j : Fin r) :
    ‖BBlockWcol (optWG (rk := alignedRk M r) β) β (alignedR M r) j‖ ^ 2 = Sagg β j := by
  rw [← real_inner_self_eq_norm_sq, inner_self_BBlockWcol_aligned β h1 j]

/-- `c_j` is an eigenvector of `A⋆` at `1 + S_j`: `A⋆ c_j = (C Cᵀ + I) c_j = (S_j + 1) c_j`,
because `Cᵀ c_j` is the `j`-th column of `Cᵀ C = diag S`. -/
theorem mulVec_BBlockWcol_aligned (β : Fin M → Fin r → ℝ) (h1 : ∀ i j, β i j ^ 2 < 1)
    (j : Fin r) :
    ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r) *ᵥ
        WithLp.ofLp (BBlockWcol (optWG (rk := alignedRk M r) β) β (alignedR M r) j)
      = (1 + Sagg β j) •
          WithLp.ofLp (BBlockWcol (optWG (rk := alignedRk M r) β) β (alignedR M r) j) := by
  set C := BBlockW (optWG (rk := alignedRk M r) β) β (alignedR M r) with hC
  have hsingle : (Pi.single j (Sagg β j) : Fin r → ℝ)
      = Sagg β j • (Pi.single j 1 : Fin r → ℝ) := by
    funext k
    by_cases h : j = k
    · rw [h]
      simp
    · simp [Ne.symm h]
  rw [ofLp_BBlockWcol, ABlockW_optWG_aligned_eq β h1, Matrix.add_mulVec, Matrix.one_mulVec,
    Matrix.mulVec_mulVec, Matrix.mul_assoc, transpose_BBlockW_optWG_aligned β h1,
    ← Matrix.mulVec_mulVec, Matrix.diagonal_mulVec_single, mul_one, hsingle,
    Matrix.mulVec_smul, add_smul, one_smul, add_comm]

/-- **Item 2.6**: at a supercritical component the projector at the sorted index `j` is the
rank-one projector on `c_j / ‖c_j‖`, so `‖P_j y‖² = ⟪c_j, y⟫² / S_j`. -/
theorem normSq_specProjIdx_ABlockW_aligned_of_sep (β : Fin M → Fin r → ℝ)
    (h0 : ∀ i j, 0 ≤ β i j) (h1 : ∀ i j, β i j < 1) (hsep : SaggSep β) (hM : 0 < M) (j : Fin r)
    (hpos : 0 < Sagg β j) (y : EuclideanSpace ℝ (Fin (rtot (alignedRk M r)))) :
    ‖specProjIdx (ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))
        (isHermitian_ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r)) (j : ℕ) y‖ ^ 2
      = ⟪BBlockWcol (optWG (rk := alignedRk M r) β) β (alignedR M r) j, y⟫_ℝ ^ 2
          / Sagg β j := by
  classical
  have hsq : ∀ i j, β i j ^ 2 < 1 := fun i j => by nlinarith [h0 i j, h1 i j]
  set A := ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r) with hA
  set hAH := isHermitian_ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r) with hAHdef
  set cj := BBlockWcol (optWG (rk := alignedRk M r) β) β (alignedR M r) j with hcj
  set S := Sagg β j with hSdef
  have hjc : (j : ℕ) < Fintype.card (Fin (rtot (alignedRk M r))) :=
    lt_card_rtot_alignedRk hM j.isLt
  have hsqrtpos : 0 < Real.sqrt S := Real.sqrt_pos.mpr hpos
  have hsqrt : Real.sqrt S ^ 2 = S := Real.sq_sqrt hpos.le
  -- the unit eigenvector at `1 + S`
  have hnormcj : ‖cj‖ = Real.sqrt S := by
    rw [← Real.sqrt_sq (norm_nonneg cj), norm_sq_BBlockWcol_aligned β hsq j]
  set u : EuclideanSpace ℝ (Fin (rtot (alignedRk M r))) := (Real.sqrt S)⁻¹ • cj with hu
  have hnormu : ‖u‖ = 1 := by
    rw [hu, norm_smul, Real.norm_eq_abs, abs_of_pos (inv_pos.mpr hsqrtpos), hnormcj,
      inv_mul_cancel₀ (ne_of_gt hsqrtpos)]
  have heigcj : toOp A cj = (1 + S) • cj := by
    have hmv : A *ᵥ WithLp.ofLp cj = (1 + S) • WithLp.ofLp cj :=
      mulVec_BBlockWcol_aligned β hsq j
    change WithLp.toLp 2 (A *ᵥ WithLp.ofLp cj) = (1 + S) • cj
    rw [hmv]
    rfl
  have heigu : toOp A u = (1 + S) • u := by
    rw [hu, map_smul, heigcj, smul_comm]
  have humem : u ∈ specSpace A (eigSetIdx A hAH (j : ℕ)) := by
    rw [eigSetIdx_ABlockW_optWG_aligned β h0 h1 hsep.1 hM j, specSpace_singleton]
    exact Module.End.mem_eigenspace_iff.mpr heigu
  have hproju : specProjIdx A hAH (j : ℕ) u = u := by
    change (specSpace A (eigSetIdx A hAH (j : ℕ))).starProjection u = u
    exact Submodule.starProjection_eq_self_iff.mpr humem
  -- the rank-one collapse
  have hsimple : SimpleSpec A hAH ((j : ℕ) + 1) :=
    simpleSpec_ABlockW_optWG_aligned_of_sep β h0 h1 hsep hM j hpos
  have hrank := fun w => specProjIdx_eq_rankOne hsimple (Nat.lt_succ_self (j : ℕ)) w
  have hvnorm : ‖vEig A hAH (j : ℕ)‖ = 1 := norm_vEig A hAH hjc
  have hut : u = ⟪vEig A hAH (j : ℕ), u⟫_ℝ • vEig A hAH (j : ℕ) := hproju.symm.trans (hrank u)
  set t := ⟪vEig A hAH (j : ℕ), u⟫_ℝ with ht
  have habs : |t| = 1 := by
    have h := hnormu
    rw [hut, norm_smul, Real.norm_eq_abs, hvnorm, mul_one] at h
    exact h
  have ht2 : t ^ 2 = 1 := by rw [← sq_abs, habs, one_pow]
  -- the two sides
  have hPy : ‖specProjIdx A hAH (j : ℕ) y‖ ^ 2 = ⟪vEig A hAH (j : ℕ), y⟫_ℝ ^ 2 := by
    rw [hrank y, norm_smul, mul_pow, Real.norm_eq_abs, sq_abs, hvnorm, one_pow, mul_one]
  have hcju : cj = Real.sqrt S • u := by
    rw [hu, smul_smul, mul_inv_cancel₀ (ne_of_gt hsqrtpos), one_smul]
  have hinner : ⟪cj, y⟫_ℝ = Real.sqrt S * (t * ⟪vEig A hAH (j : ℕ), y⟫_ℝ) := by
    rw [hcju, real_inner_smul_left, hut, real_inner_smul_left]
  rw [hPy, hinner, mul_pow, mul_pow, hsqrt, ht2, one_mul, mul_comm, mul_div_assoc,
    div_self (ne_of_gt hpos), mul_one]

/-- **Item 2.6** at the paper's hypothesis (b). -/
theorem normSq_specProjIdx_ABlockW_aligned (β : Fin M → Fin r → ℝ) (h0 : ∀ i j, 0 ≤ β i j)
    (h1 : ∀ i j, β i j < 1) (hS : StrictAnti (Sagg β)) (hM : 0 < M) (j : Fin r)
    (hpos : 0 < Sagg β j) (y : EuclideanSpace ℝ (Fin (rtot (alignedRk M r)))) :
    ‖specProjIdx (ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))
        (isHermitian_ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r)) (j : ℕ) y‖ ^ 2
      = ⟪BBlockWcol (optWG (rk := alignedRk M r) β) β (alignedR M r) j, y⟫_ℝ ^ 2
          / Sagg β j :=
  normSq_specProjIdx_ABlockW_aligned_of_sep β h0 h1 (Sagg_sep_of_strictAnti hS) hM j hpos y

/-- **Item 2.7**: the value the component clause converges to. At `y = c_j` the quotient
`‖P_j c_j‖² / (1 + S_j)` is `S_j / (S_j + 1)`, because `‖c_j‖² = S_j`. -/
theorem compLimit_aligned_of_sep (β : Fin M → Fin r → ℝ) (h0 : ∀ i j, 0 ≤ β i j)
    (h1 : ∀ i j, β i j < 1) (hsep : SaggSep β) (hM : 0 < M) (j : Fin r)
    (hpos : 0 < Sagg β j) :
    ‖specProjIdx (ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))
          (isHermitian_ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r)) (j : ℕ)
          (BBlockWcol (optWG (rk := alignedRk M r) β) β (alignedR M r) j)‖ ^ 2
        / (1 + Sagg β j)
      = Sagg β j / (Sagg β j + 1) := by
  have hsq : ∀ i j, β i j ^ 2 < 1 := fun i j => by nlinarith [h0 i j, h1 i j]
  rw [normSq_specProjIdx_ABlockW_aligned_of_sep β h0 h1 hsep hM j hpos,
    inner_self_BBlockWcol_aligned β hsq j, sq, mul_div_assoc,
    div_self (ne_of_gt hpos), mul_one, add_comm]

/-- **Item 2.7** at the paper's hypothesis (b). -/
theorem compLimit_aligned (β : Fin M → Fin r → ℝ) (h0 : ∀ i j, 0 ≤ β i j)
    (h1 : ∀ i j, β i j < 1) (hS : StrictAnti (Sagg β)) (hM : 0 < M) (j : Fin r)
    (hpos : 0 < Sagg β j) :
    ‖specProjIdx (ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))
          (isHermitian_ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r)) (j : ℕ)
          (BBlockWcol (optWG (rk := alignedRk M r) β) β (alignedR M r) j)‖ ^ 2
        / (1 + Sagg β j)
      = Sagg β j / (Sagg β j + 1) :=
  compLimit_aligned_of_sep β h0 h1 (Sagg_sep_of_strictAnti hS) hM j hpos

end Projector

/-! ### 4. `S` is antitone inside the model (item 6.3)

Finding F1 of `notes/archive/rankr_D4_plan.md`: the paper's hypothesis (b) `S_j ≠ S_k` follows from
the sorted form of hypothesis (a) at every index that matters (not from the unordered (a) of
the paper; see the header). Inside the Lean model the field
`SpikedModelR.hθanti` orders every `θ_i` strictly, and
`β²/(1-β²) = svdTerm θ c = (θ⁴-c)/(θ²+c)` above the threshold is monotone in `θ ≥ 0`. So
`Sagg β` is antitone, and strictly decreasing from any index that carries `S_j > 0`. -/

section Monotone

/-- `svdTerm` is monotone in `θ` on `[0, ∞)`. Above the threshold the cross-multiplied
inequality is `a²b²(b²-a²) + c(b⁴-a⁴) + c(b²-a²) ≥ 0`; below it the value is `0`. -/
theorem svdTerm_mono {a b c : ℝ} (hc : 0 < c) (ha : 0 ≤ a) (hab : a ≤ b) :
    Scalars.svdTerm a c ≤ Scalars.svdTerm b c := by
  have hb : 0 ≤ b := le_trans ha hab
  unfold Scalars.svdTerm
  split_ifs with hA hB hB
  · have hda : (0 : ℝ) < a ^ 2 + c := by nlinarith [sq_nonneg a]
    have hdb : (0 : ℝ) < b ^ 2 + c := by nlinarith [sq_nonneg b]
    rw [div_le_div_iff₀ hda hdb]
    have h2 : a ^ 2 ≤ b ^ 2 := by nlinarith
    have h4 : a ^ 4 ≤ b ^ 4 := by
      nlinarith [mul_le_mul h2 h2 (sq_nonneg a) (sq_nonneg b)]
    have hab2 : (0 : ℝ) ≤ a ^ 2 * b ^ 2 := by positivity
    nlinarith [mul_nonneg hab2 (sub_nonneg.mpr h2), mul_nonneg hc.le (sub_nonneg.mpr h4),
      mul_nonneg hc.le (sub_nonneg.mpr h2)]
  · exfalso
    have h2 : a ^ 2 ≤ b ^ 2 := by nlinarith
    have h4 : a ^ 4 ≤ b ^ 4 := by nlinarith [sq_nonneg a, sq_nonneg b]
    exact hB (lt_of_lt_of_le hA h4)
  · have hdb : (0 : ℝ) < b ^ 2 + c := by nlinarith [sq_nonneg b]
    exact div_nonneg (by linarith) hdb.le
  · exact le_rfl

/-- `svdTerm` is strictly monotone in `θ` wherever it is positive: below a supercritical `b`
every `a < b` gives a strictly smaller value. -/
theorem svdTerm_lt_of_pos {a b c : ℝ} (hc : 0 < c) (ha : 0 ≤ a) (hab : a < b)
    (hpos : 0 < Scalars.svdTerm b c) : Scalars.svdTerm a c < Scalars.svdTerm b c := by
  have hb : 0 ≤ b := le_trans ha hab.le
  have hB : c < b ^ 4 := by
    by_contra hcon
    rw [Scalars.svdTerm, if_neg hcon] at hpos
    exact lt_irrefl 0 hpos
  have hdb : (0 : ℝ) < b ^ 2 + c := by nlinarith [sq_nonneg b]
  unfold Scalars.svdTerm at hpos ⊢
  rw [if_pos hB] at hpos ⊢
  by_cases hA : c < a ^ 4
  · rw [if_pos hA]
    have hda : (0 : ℝ) < a ^ 2 + c := by nlinarith [sq_nonneg a]
    rw [div_lt_div_iff₀ hda hdb]
    have h2 : a ^ 2 < b ^ 2 := by nlinarith
    have hsum : (0 : ℝ) < b ^ 2 + a ^ 2 := by nlinarith [sq_nonneg a]
    have h4 : a ^ 4 < b ^ 4 := by
      nlinarith [mul_pos (sub_pos.mpr h2) hsum]
    have hab2 : (0 : ℝ) ≤ a ^ 2 * b ^ 2 := by positivity
    nlinarith [mul_nonneg hab2 (sub_nonneg.mpr h2.le), mul_pos hc (sub_pos.mpr h4),
      mul_pos hc (sub_pos.mpr h2)]
  · rw [if_neg hA]
    exact hpos

/-- `S_j` in the paper's data `(θ, c)`: `S_j = ∑_i svdTerm θ_ij c_i` (`main_paper.tex:2395`).
The rank-`r` twin of `Sval_beta` (`Scalars.lean`). -/
theorem Sagg_beta (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i) (j : Fin r) :
    Sagg (fun i j => beta (θ i j) (c i)) j = ∑ i, Scalars.svdTerm (θ i j) (c i) := by
  rw [Sagg]
  exact Finset.sum_congr rfl fun i _ => by
    rw [Scalars.beta_sq, Scalars.betaSq_div_one_sub (hc i)]

/-- **Item 6.3**: `Sagg β` is antitone once every table orders its spikes. -/
theorem Sagg_antitone {θ : Fin M → Fin r → ℝ} {c : Fin M → ℝ} {β : Fin M → Fin r → ℝ}
    (hc : ∀ i, 0 < c i) (hθpos : ∀ i j, 0 ≤ θ i j) (hθanti : ∀ i, Antitone (θ i))
    (hβdef : ∀ i j, β i j = beta (θ i j) (c i)) : Antitone (Sagg β) := by
  have hβ : β = fun i j => beta (θ i j) (c i) := funext fun i => funext fun j => hβdef i j
  subst hβ
  intro j k hjk
  rw [Sagg_beta θ c hc, Sagg_beta θ c hc]
  exact Finset.sum_le_sum fun i _ =>
    svdTerm_mono (hc i) (hθpos i k) (hθanti i hjk)

/-- **Item 6.3**: `Sagg β` decreases strictly from any component that carries `S_j > 0`. With
`Sagg_antitone` this is the paper's hypothesis (b) `S_j ≠ S_k`, so `hS : StrictAnti (Sagg β)`
is not an extra assumption at the indices where the component clause has content. -/
theorem Sagg_lt_of_pos {θ : Fin M → Fin r → ℝ} {c : Fin M → ℝ} {β : Fin M → Fin r → ℝ}
    (hc : ∀ i, 0 < c i) (hθpos : ∀ i j, 0 ≤ θ i j) (hθanti : ∀ i, StrictAnti (θ i))
    (hβdef : ∀ i j, β i j = beta (θ i j) (c i)) {j k : Fin r} (hjk : j < k)
    (hpos : 0 < Sagg β j) : Sagg β k < Sagg β j := by
  have hβ : β = fun i j => beta (θ i j) (c i) := funext fun i => funext fun j => hβdef i j
  subst hβ
  rw [Sagg_beta θ c hc] at hpos
  rw [Sagg_beta θ c hc, Sagg_beta θ c hc]
  have hex : ∃ i : Fin M, 0 < Scalars.svdTerm (θ i j) (c i) := by
    by_contra hcon
    have hle : ∑ i, Scalars.svdTerm (θ i j) (c i) ≤ 0 :=
      Finset.sum_nonpos fun i _ => not_lt.mp fun h => hcon ⟨i, h⟩
    linarith
  obtain ⟨i₀, hi₀⟩ := hex
  refine Finset.sum_lt_sum
    (fun i _ => svdTerm_mono (hc i) (hθpos i k) (le_of_lt (hθanti i hjk)))
    ⟨i₀, Finset.mem_univ i₀, svdTerm_lt_of_pos (hc i₀) (hθpos i₀ k) (hθanti i₀ hjk) hi₀⟩

end Monotone

namespace UnalignedModelR

open MeasureTheory

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-- `Sagg β` is antitone inside the model: `SpikedModelR.hθanti` and `SpikedModelR.hθnn` supply
the two hypotheses of `Sagg_antitone`. -/
theorem Sagg_antitone_of_model (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) : Antitone (Sagg β) :=
  Sagg_antitone hc (fun i j => (m.tbl i).hθnn j)
    (fun i => (m.tbl i).hθanti.antitone) hβdef

/-- `Sagg β` decreases strictly from any supercritical component, inside the model. -/
theorem Sagg_lt_of_pos_of_model (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) {j k : Fin r} (hjk : j < k)
    (hpos : 0 < Sagg β j) : Sagg β k < Sagg β j :=
  Sagg_lt_of_pos hc (fun i j => (m.tbl i).hθnn j) (fun i => (m.tbl i).hθanti) hβdef
    hjk hpos

/-- **The `hS`-free separation inside the model.** `SaggSep β` is everything `padSagg_lt`
reads, so the gaps and the simplicity of item 2.5 are available with no hypothesis (b): the
model's own spike order supplies them. -/
theorem saggSep_of_model (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) : SaggSep β :=
  ⟨m.Sagg_antitone_of_model c β hc hβdef,
    fun _ _ hab hp => m.Sagg_lt_of_pos_of_model c β hc hβdef hab hp⟩

end UnalignedModelR

end StackedSVD
