/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.Duality
import StackedSVD.LinAlg.SpecIdxPerturb
import StackedSVD.LinAlg.Eigen

/-!
# Duality at a sorted index between the two Gram matrices

Stage E2 of `notes/archive/rankr_TrackE_plan.md` (section 2.3, step 6). The column split of the
heteroscedastic rank-`r` model puts the rank-`r` update on the `n` side, that is on `X Xᵀ`,
while `overlapIdx` reads the spectral projector of the `d` side `Xᵀ X` at one sorted index.
Four lemmas cross the gap, with no simplicity hypothesis: the spectral projector at the index
replaces the single eigenvector of the rank-1 mirror `RMT/Het/Duality.lean`.

* `eigVal_gram_comm`: `λ_k(Xᵀ X) = λ_k(X Xᵀ)` for every `k < min p q`. Mirror:
  `Het.lamMax_gram_comm`, which is the case `k = 0`. Route: `Matrix.charpoly_mul_comm_of_le`
  gives `(X Xᵀ).charpoly = X^{p-q} (Xᵀ X).charpoly` when `q ≤ p`, so the sorted spectrum of
  `X Xᵀ` is that of `Xᵀ X` padded by zeros (`eigenvalues₀_eq_of_charpoly`,
  `antitone_padSpec`, `prod_X_sub_C_padSpec` of `LinAlg/Eigen.lean`); the case `p ≤ q` is the
  same at `Xᵀ`.
* `specProjIdx_mul_transpose_mulVec`: `X` intertwines the two projectors, `P'_ℓ X = X P_ℓ`.
* `normSq_mulVec_specProjIdx`: `‖X P_ℓ y‖² = λ_ℓ ‖P_ℓ y‖²`.
* `overlapIdx_eq_normSq_specProjIdx_div`: **duality at a sorted index**,
  `overlapIdx X ℓ y = ‖P'_ℓ (X y)‖² / λ_ℓ(X Xᵀ)` when `λ_ℓ > 0`. Mirror:
  `Het.topProj_transpose_eq`. Its right side is `compFun ℓ k` of `LinAlg/SpecIdxPerturb.lean`
  at `C = X Xᵀ` and `Y = X V`.
* `simpleIdx_gram_comm`: simplicity at a positive sorted index transfers between the two Gram
  matrices. Mirror: `topSimple_transpose_mul_iff` (`LinAlg/TopProjPerturb.lean`).

Paper: `main_paper.tex` lines 1409 to 1440 (rank 1) and 2294 to 2331 (rank `r`).

Everything here is model-free: no `UnalignedModelR` appears. Stage E1 lands
`stackXW_mulVec_colVecG`, and stage E7 composes the two.
-/

open scoped Matrix InnerProductSpace

namespace StackedSVD
namespace Het

variable {p q : ℕ}

/-! ### 1. Small helpers -/

/-- `(X y) ⬝ᵥ z = y ⬝ᵥ (Xᵀ z)`, the adjoint identity in `dotProduct` form. A copy of the
private helper of `RMT/Het/Duality.lean`. -/
private theorem dotProduct_mulVec_transpose' (X : Matrix (Fin p) (Fin q) ℝ)
    (y : Fin q → ℝ) (z : Fin p → ℝ) : (X *ᵥ y) ⬝ᵥ z = y ⬝ᵥ (Xᵀ *ᵥ z) := by
  simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply, Finset.sum_mul, Finset.mul_sum]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun k _ => Finset.sum_congr rfl fun i _ => by ring

/-- `‖x‖²` as a dot product. A copy of the private helper of `RMT/Het/Duality.lean`. -/
private theorem dotProduct_self_eq_norm_sq' {m : ℕ} (x : EuclideanSpace ℝ (Fin m)) :
    WithLp.ofLp x ⬝ᵥ WithLp.ofLp x = ‖x‖ ^ 2 := by
  have h := real_inner_eq_dotProduct x x
  rw [real_inner_self_eq_norm_sq] at h
  exact h.symm

/-- The sorted eigenvalues of `Xᵀ X` are nonnegative. -/
private theorem eigenvalues₀_transpose_mul_self_nonneg (X : Matrix (Fin p) (Fin q) ℝ)
    (k : Fin (Fintype.card (Fin q))) :
    0 ≤ (isHermitian_transpose_mul_self X).eigenvalues₀ k := by
  have hpsd : (Xᵀ * X).PosSemidef := by
    simpa using Matrix.posSemidef_conjTranspose_mul_self X
  rw [← eigenvalues_eigIdx]
  exact hpsd.eigenvalues_nonneg _

/-- The sorted eigenvalues of `X Xᵀ` are nonnegative. -/
private theorem eigenvalues₀_mul_transpose_self_nonneg (X : Matrix (Fin p) (Fin q) ℝ)
    (k : Fin (Fintype.card (Fin p))) :
    0 ≤ (isHermitian_mul_transpose_self X).eigenvalues₀ k := by
  have hpsd : (X * Xᵀ).PosSemidef := by
    simpa using Matrix.posSemidef_self_mul_conjTranspose X
  rw [← eigenvalues_eigIdx]
  exact hpsd.eigenvalues_nonneg _

/-- The characteristic polynomial of a Hermitian matrix is the product over its sorted
eigenvalues. `charpoly_eq_prod_of_conj` (`LinAlg/Eigen.lean`) reindexed by `eigIdx`. -/
private theorem charpoly_eq_prod_eigenvalues₀ {n : ℕ} {A : Matrix (Fin n) (Fin n) ℝ}
    (hA : A.IsHermitian) :
    A.charpoly = ∏ k : Fin (Fintype.card (Fin n)),
      (Polynomial.X - Polynomial.C (hA.eigenvalues₀ k)) := by
  rw [charpoly_eq_prod_of_conj (Unitary.coe_star_mul_self _) (spectral_conj hA)]
  exact (Fintype.prod_equiv (eigIdx n) _ _ fun k => by rw [eigenvalues_eigIdx]).symm

/-! ### 2. The sorted eigenvalues of the two Gram matrices agree -/

/-- The case `q ≤ p` of `eigVal_gram_comm`: the sorted spectrum of `X Xᵀ` is that of `Xᵀ X`
padded by `p - q` zeros. The template is
`eigenvalues₀_mul_transpose_of_transpose_mul_diagonal` (`LinAlg/Eigen.lean`). -/
private theorem eigenvalues₀_mul_transpose_self_eq_pad (X : Matrix (Fin p) (Fin q) ℝ)
    (hqp : q ≤ p) (k : Fin (Fintype.card (Fin p))) :
    (isHermitian_mul_transpose_self X).eigenvalues₀ k
      = if h : (k : ℕ) < Fintype.card (Fin q)
          then (isHermitian_transpose_mul_self X).eigenvalues₀ ⟨k, h⟩ else 0 := by
  classical
  have hle : Fintype.card (Fin q) ≤ Fintype.card (Fin p) := by simpa using hqp
  have hchar : (X * Xᵀ).charpoly
      = ∏ l : Fin (Fintype.card (Fin p)), (Polynomial.X - Polynomial.C
          (if h : (l : ℕ) < Fintype.card (Fin q)
            then (isHermitian_transpose_mul_self X).eigenvalues₀ ⟨l, h⟩ else 0)) := by
    rw [prod_X_sub_C_padSpec hle, Matrix.charpoly_mul_comm_of_le X Xᵀ hle,
      charpoly_eq_prod_eigenvalues₀ (isHermitian_transpose_mul_self X)]
  exact congrFun (eigenvalues₀_eq_of_charpoly (isHermitian_mul_transpose_self X)
    (antitone_padSpec (isHermitian_transpose_mul_self X).eigenvalues₀_antitone
      (eigenvalues₀_transpose_mul_self_nonneg X)) hchar) k

/-- **The sorted eigenvalues of `Xᵀ X` and `X Xᵀ` agree at every index below `min p q`.**
Mirror: `Het.lamMax_gram_comm` (`RMT/Het/Duality.lean`), which is the case `k = 0`. Beyond
`min p q` the longer list carries zeros, so the statement is false there in general. -/
theorem eigVal_gram_comm (X : Matrix (Fin p) (Fin q) ℝ) {k : ℕ} (hk : k < min p q) :
    eigVal (Xᵀ * X) (isHermitian_transpose_mul_self X) k
      = eigVal (X * Xᵀ) (isHermitian_mul_transpose_self X) k := by
  have hkq : k < Fintype.card (Fin q) := by rw [Fintype.card_fin]; omega
  have hkp : k < Fintype.card (Fin p) := by rw [Fintype.card_fin]; omega
  rw [eigVal_eq _ _ hkq, eigVal_eq _ _ hkp]
  rcases le_total q p with hqp | hpq
  · rw [eigenvalues₀_mul_transpose_self_eq_pad X hqp ⟨k, hkp⟩, dif_pos hkq]
  · have h := eigenvalues₀_mul_transpose_self_eq_pad Xᵀ hpq ⟨k, hkq⟩
    rw [dif_pos hkp] at h
    have e1 : (isHermitian_mul_transpose_self Xᵀ).eigenvalues₀
        = (isHermitian_transpose_mul_self X).eigenvalues₀ :=
      eigenvalues₀_congr_mat _ _ (by rw [Matrix.transpose_transpose])
    have e2 : (isHermitian_transpose_mul_self Xᵀ).eigenvalues₀
        = (isHermitian_mul_transpose_self X).eigenvalues₀ :=
      eigenvalues₀_congr_mat _ _ (by rw [Matrix.transpose_transpose])
    rw [e1, e2] at h
    exact h

/-! ### 3. `X` intertwines the two spectral projectors -/

/-- **`P'_λ X = X P_λ`.** For `ℓ < min p q`, the spectral projector of `X Xᵀ` at the sorted
index `ℓ` applied to `X y` is `X` applied to the projector of `Xᵀ X` at the same index.
Route: `X P y` lies in the `λ_ℓ`-eigenspace of `X Xᵀ`, and `X (1 - P) y` is orthogonal to
it, which characterizes the orthogonal projection
(`Submodule.eq_starProjection_of_mem_of_inner_eq_zero`). -/
theorem specProjIdx_mul_transpose_mulVec (X : Matrix (Fin p) (Fin q) ℝ) {ℓ : ℕ}
    (hℓ : ℓ < min p q) (y : EuclideanSpace ℝ (Fin q)) :
    specProjIdx (X * Xᵀ) (isHermitian_mul_transpose_self X) ℓ
        (WithLp.toLp 2 (X *ᵥ WithLp.ofLp y))
      = WithLp.toLp 2 (X *ᵥ WithLp.ofLp
          (specProjIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) ℓ y)) := by
  have hkq : ℓ < Fintype.card (Fin q) := by rw [Fintype.card_fin]; omega
  have hkp : ℓ < Fintype.card (Fin p) := by rw [Fintype.card_fin]; omega
  set lam := (isHermitian_transpose_mul_self X).eigenvalues₀ ⟨ℓ, hkq⟩ with hlam
  have hlamB : (isHermitian_mul_transpose_self X).eigenvalues₀ ⟨ℓ, hkp⟩ = lam := by
    have h := eigVal_gram_comm X hℓ
    rw [eigVal_eq _ _ hkq, eigVal_eq _ _ hkp] at h
    exact h.symm
  -- the two spectral subspaces are the eigenspaces at `lam`
  have hKA : specSpace (Xᵀ * X) (eigSetIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) ℓ)
      = Module.End.eigenspace (toOp (Xᵀ * X)) lam := by
    rw [eigSetIdx_eq_singleton _ _ hkq, specSpace_singleton]
  have hKB : specSpace (X * Xᵀ) (eigSetIdx (X * Xᵀ) (isHermitian_mul_transpose_self X) ℓ)
      = Module.End.eigenspace (toOp (X * Xᵀ)) lam := by
    rw [eigSetIdx_eq_singleton _ _ hkp, hlamB, specSpace_singleton]
  set P := specProjIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) ℓ with hP
  -- `P y` is in the eigenspace of `Xᵀ X`
  have hPy : (Xᵀ * X) *ᵥ WithLp.ofLp (P y) = lam • WithLp.ofLp (P y) := by
    have hmem : P y
        ∈ specSpace (Xᵀ * X) (eigSetIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) ℓ) :=
      Submodule.starProjection_apply_mem _ y
    rw [hKA, R4.mem_eigenspace_iff'] at hmem
    exact hmem
  change Submodule.starProjection
    (specSpace (X * Xᵀ) (eigSetIdx (X * Xᵀ) (isHermitian_mul_transpose_self X) ℓ)) _ = _
  refine Submodule.eq_starProjection_of_mem_of_inner_eq_zero ?_ ?_
  · -- `X (P y)` is in the eigenspace of `X Xᵀ`
    rw [hKB, R4.mem_eigenspace_iff']
    change (X * Xᵀ) *ᵥ (X *ᵥ WithLp.ofLp (P y)) = lam • (X *ᵥ WithLp.ofLp (P y))
    rw [Matrix.mulVec_mulVec, Matrix.mul_assoc, ← Matrix.mulVec_mulVec, hPy,
      Matrix.mulVec_smul]
  · intro w hw
    rw [hKB, R4.mem_eigenspace_iff'] at hw
    -- `Xᵀ w` is in the eigenspace of `Xᵀ X`
    have hXw : (Xᵀ * X) *ᵥ (Xᵀ *ᵥ WithLp.ofLp w) = lam • (Xᵀ *ᵥ WithLp.ofLp w) := by
      rw [Matrix.mulVec_mulVec, Matrix.mul_assoc, ← Matrix.mulVec_mulVec, hw,
        Matrix.mulVec_smul]
    have hXwmem : (WithLp.toLp 2 (Xᵀ *ᵥ WithLp.ofLp w) : EuclideanSpace ℝ (Fin q))
        ∈ specSpace (Xᵀ * X) (eigSetIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) ℓ) := by
      rw [hKA, R4.mem_eigenspace_iff']
      exact hXw
    have horth : ⟪y - P y,
        (WithLp.toLp 2 (Xᵀ *ᵥ WithLp.ofLp w) : EuclideanSpace ℝ (Fin q))⟫_ℝ = 0 :=
      Submodule.starProjection_inner_eq_zero y _ hXwmem
    rw [real_inner_eq_dotProduct] at horth ⊢
    change (X *ᵥ WithLp.ofLp y - X *ᵥ WithLp.ofLp (P y)) ⬝ᵥ WithLp.ofLp w = 0
    rw [← Matrix.mulVec_sub, dotProduct_mulVec_transpose']
    exact horth

/-! ### 4. The norm of `X` on the eigenspace -/

/-- **`‖X z‖² = λ ‖z‖²` on the `λ`-eigenspace of `Xᵀ X`**, read at the projector: for
`z = P_ℓ y`, `‖X z‖² = ⟪z, Xᵀ X z⟫ = λ_ℓ ‖z‖²`. Out of range both sides are `0`, so no
bound on `ℓ` is needed. -/
theorem normSq_mulVec_specProjIdx (X : Matrix (Fin p) (Fin q) ℝ) (ℓ : ℕ)
    (y : EuclideanSpace ℝ (Fin q)) :
    ‖(WithLp.toLp 2 (X *ᵥ WithLp.ofLp
        (specProjIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) ℓ y))
        : EuclideanSpace ℝ (Fin p))‖ ^ 2
      = eigVal (Xᵀ * X) (isHermitian_transpose_mul_self X) ℓ
          * ‖specProjIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) ℓ y‖ ^ 2 := by
  set P := specProjIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) ℓ with hP
  by_cases hkq : ℓ < Fintype.card (Fin q)
  · set lam := (isHermitian_transpose_mul_self X).eigenvalues₀ ⟨ℓ, hkq⟩ with hlam
    have hPy : (Xᵀ * X) *ᵥ WithLp.ofLp (P y) = lam • WithLp.ofLp (P y) := by
      have hmem : P y
          ∈ specSpace (Xᵀ * X) (eigSetIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) ℓ) :=
        Submodule.starProjection_apply_mem _ y
      rw [eigSetIdx_eq_singleton _ _ hkq, specSpace_singleton, R4.mem_eigenspace_iff'] at hmem
      exact hmem
    rw [eigVal_eq _ _ hkq, ← dotProduct_self_eq_norm_sq', ← dotProduct_self_eq_norm_sq']
    change (X *ᵥ WithLp.ofLp (P y)) ⬝ᵥ (X *ᵥ WithLp.ofLp (P y))
      = lam * (WithLp.ofLp (P y) ⬝ᵥ WithLp.ofLp (P y))
    rw [dotProduct_mulVec_transpose', Matrix.mulVec_mulVec, hPy, dotProduct_smul, smul_eq_mul]
  · rw [hP, specProjIdx_apply_of_not_lt _ _ hkq, eigVal_of_le _ _ (not_lt.mp hkq)]
    simp

/-! ### 5. Duality at a sorted index -/

/-- A positive sorted eigenvalue of `Xᵀ X` sits below `min p q`: at or beyond `q` the value is
the junk `0`, and at or beyond `p` it is `0` because the rank is at most `p`. So `hℓ` follows
from `hpos` in `overlapIdx_eq_normSq_specProjIdx_div` and `simpleIdx_gram_comm`; both call this
lemma, so neither takes `hℓ` as a hypothesis. -/
theorem lt_min_of_eigVal_transpose_mul_self_pos (X : Matrix (Fin p) (Fin q) ℝ) {ℓ : ℕ}
    (hpos : 0 < eigVal (Xᵀ * X) (isHermitian_transpose_mul_self X) ℓ) : ℓ < min p q := by
  by_contra hcon
  rw [not_lt] at hcon
  rcases le_or_gt (Fintype.card (Fin q)) ℓ with hq | hq
  · rw [eigVal_of_le _ _ hq] at hpos
    exact lt_irrefl _ hpos
  · rw [eigVal_eq _ _ hq] at hpos
    have hq' : ℓ < q := by simpa using hq
    have hp : p ≤ ℓ := by omega
    have hrank : (Xᵀ * X).rank ≤ p := by
      rw [Matrix.rank_transpose_mul_self]; exact Matrix.rank_le_height X
    have hzero := eigenvalues₀_eq_zero_of_rank_le (isHermitian_transpose_mul_self X)
      (eigenvalues₀_transpose_mul_self_nonneg X) hrank (l := ⟨ℓ, hq⟩) hp
    rw [hzero] at hpos
    exact lt_irrefl _ hpos

/-- **Duality at a sorted index** (plan step 6; rank-1 mirror `Het.topProj_transpose_eq`).
The `d`-side overlap at the sorted index `ℓ` is the `n`-side projector norm of `X y`, divided
by the shared eigenvalue: `overlapIdx X ℓ y = ‖P'_ℓ (X y)‖² / λ_ℓ(X Xᵀ)`. `hpos` gives the
range bound `ℓ < min p q` (`lt_min_of_eigVal_transpose_mul_self_pos`), which
`eigVal_gram_comm` and `specProjIdx_mul_transpose_mulVec` need, and it also makes the
division legal. -/
theorem overlapIdx_eq_normSq_specProjIdx_div (X : Matrix (Fin p) (Fin q) ℝ) {ℓ : ℕ}
    (hpos : 0 < eigVal (Xᵀ * X) (isHermitian_transpose_mul_self X) ℓ)
    (y : EuclideanSpace ℝ (Fin q)) :
    overlapIdx X ℓ y
      = ‖specProjIdx (X * Xᵀ) (isHermitian_mul_transpose_self X) ℓ
            (WithLp.toLp 2 (X *ᵥ WithLp.ofLp y))‖ ^ 2
        / eigVal (X * Xᵀ) (isHermitian_mul_transpose_self X) ℓ := by
  have hℓ : ℓ < min p q := lt_min_of_eigVal_transpose_mul_self_pos X hpos
  have hlam := eigVal_gram_comm X hℓ
  have hne : eigVal (X * Xᵀ) (isHermitian_mul_transpose_self X) ℓ ≠ 0 := by
    rw [← hlam]; exact ne_of_gt hpos
  rw [specProjIdx_mul_transpose_mulVec X hℓ y, normSq_mulVec_specProjIdx X ℓ y, hlam,
    overlapIdx, mul_div_cancel_left₀ _ hne]

/-! ### 6. Simplicity transfers -/

/-- `SimpleIdx` transports along a matrix equality. -/
private theorem simpleIdx_congr_mat {n : ℕ} {A B : Matrix (Fin n) (Fin n) ℝ}
    (hA : A.IsHermitian) (hB : B.IsHermitian) (h : A = B) (j : ℕ) :
    SimpleIdx A hA j ↔ SimpleIdx B hB j := by
  subst h; exact Iff.rfl

/-- One direction of `simpleIdx_gram_comm`: from `Xᵀ X` to `X Xᵀ`. -/
private theorem simpleIdx_mul_transpose_of (X : Matrix (Fin p) (Fin q) ℝ) {ℓ : ℕ}
    (hℓ : ℓ < min p q) (hpos : 0 < eigVal (Xᵀ * X) (isHermitian_transpose_mul_self X) ℓ)
    (hs : SimpleIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) ℓ) :
    SimpleIdx (X * Xᵀ) (isHermitian_mul_transpose_self X) ℓ := by
  intro k l hk hkl hcon
  have hkq : ℓ < Fintype.card (Fin q) := by rw [Fintype.card_fin]; omega
  have hkp : ℓ < Fintype.card (Fin p) := by rw [Fintype.card_fin]; omega
  have hkeq : k = ⟨ℓ, hkp⟩ := Fin.ext hk
  subst hkeq
  have hlamB : (isHermitian_mul_transpose_self X).eigenvalues₀ ⟨ℓ, hkp⟩
      = (isHermitian_transpose_mul_self X).eigenvalues₀ ⟨ℓ, hkq⟩ := by
    have h := eigVal_gram_comm X hℓ
    rw [eigVal_eq _ _ hkq, eigVal_eq _ _ hkp] at h
    exact h.symm
  have hlp : (l : ℕ) < p := by simpa using l.isLt
  by_cases hl : (l : ℕ) < q
  · -- `l` is below `min p q`: the eigenvalues of the two sides agree there
    have hlq : (l : ℕ) < Fintype.card (Fin q) := by simpa using hl
    have hlB : (isHermitian_mul_transpose_self X).eigenvalues₀ l
        = (isHermitian_transpose_mul_self X).eigenvalues₀ ⟨l, hlq⟩ := by
      have h := eigVal_gram_comm X (lt_min hlp hl)
      rw [eigVal_eq _ _ hlq, eigVal_eq _ _ l.isLt] at h
      exact h.symm
    have hne : (⟨ℓ, hkq⟩ : Fin (Fintype.card (Fin q))) ≠ ⟨l, hlq⟩ := by
      intro h
      exact hkl (Fin.ext (by simpa using congrArg Fin.val h))
    refine hs ⟨ℓ, hkq⟩ ⟨l, hlq⟩ rfl hne ?_
    rw [← hlB, ← hlamB]
    exact hcon
  · -- `l` is at or beyond `q`, so its eigenvalue on the `X Xᵀ` side is `0`
    have hrank : (X * Xᵀ).rank ≤ q := by
      rw [Matrix.rank_self_mul_transpose]; exact Matrix.rank_le_width X
    have hzero : (isHermitian_mul_transpose_self X).eigenvalues₀ l = 0 :=
      eigenvalues₀_eq_zero_of_rank_le _ (eigenvalues₀_mul_transpose_self_nonneg X) hrank
        (not_lt.mp hl)
    rw [hzero, hlamB] at hcon
    rw [eigVal_eq _ _ hkq, ← hcon] at hpos
    exact lt_irrefl _ hpos

/-- **Simplicity at a positive sorted index transfers between the two Gram matrices.**
Mirror: `topSimple_transpose_mul_iff` (`LinAlg/TopProjPerturb.lean`). Indices below
`min p q` agree by `eigVal_gram_comm`; an index at or beyond `min p q` on either side carries
the eigenvalue `0` (`eigenvalues₀_eq_zero_of_rank_le`, the rank is at most `min p q`), which
`hpos` separates from `λ_ℓ`. `hpos` also gives the range bound `ℓ < min p q`
(`lt_min_of_eigVal_transpose_mul_self_pos`). -/
theorem simpleIdx_gram_comm (X : Matrix (Fin p) (Fin q) ℝ) {ℓ : ℕ}
    (hpos : 0 < eigVal (Xᵀ * X) (isHermitian_transpose_mul_self X) ℓ) :
    SimpleIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) ℓ
      ↔ SimpleIdx (X * Xᵀ) (isHermitian_mul_transpose_self X) ℓ := by
  have hℓ : ℓ < min p q := lt_min_of_eigVal_transpose_mul_self_pos X hpos
  refine ⟨simpleIdx_mul_transpose_of X hℓ hpos, fun hs => ?_⟩
  have hℓ' : ℓ < min q p := by rw [min_comm]; exact hℓ
  have hpos' : 0 < eigVal (Xᵀᵀ * Xᵀ) (isHermitian_transpose_mul_self Xᵀ) ℓ := by
    rw [eigVal_congr_mat (isHermitian_transpose_mul_self Xᵀ) (isHermitian_mul_transpose_self X)
      (by rw [Matrix.transpose_transpose]) ℓ, ← eigVal_gram_comm X hℓ]
    exact hpos
  have hs' : SimpleIdx (Xᵀᵀ * Xᵀ) (isHermitian_transpose_mul_self Xᵀ) ℓ :=
    (simpleIdx_congr_mat (isHermitian_transpose_mul_self Xᵀ) (isHermitian_mul_transpose_self X)
      (by rw [Matrix.transpose_transpose]) ℓ).mpr hs
  exact (simpleIdx_congr_mat (isHermitian_mul_transpose_self Xᵀ)
    (isHermitian_transpose_mul_self X) (by rw [Matrix.transpose_transpose]) ℓ).mp
    (simpleIdx_mul_transpose_of Xᵀ hℓ' hpos' hs')

end Het
end StackedSVD
