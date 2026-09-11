/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.LinAlg.Eigen
import StackedSVD.RMT.R4

/-!
# The frame lemma: an approximate eigenframe against the threshold projector

STATUS 2026-08-30: written for task U3 of `notes/archive/plan_subspacelaw.md` (Route P, section 1.3
items 3 to 5, section 2 "New", section 4 risk 1). Model free and deterministic: one real
symmetric matrix `S : Matrix (Fin p) (Fin p) ℝ`, a threshold `τ`, and a family of `s`
approximate eigenvectors.

## Content

1. Elementary vector inequalities for entrywise bounds (`abs_dotProduct_mulVec_le`,
   `sum_sq_mulVec_le`), the injectivity of a matrix whose Gram is close to `1`
   (`mulVec_eq_zero_of_gram_close`, `linearIndependent_of_gram_close`), and the square
   frame lemma `sq_frame_close` : if `B : Matrix ι (Fin n) ℝ`, `card ι = n`, and
   `|BᵀB - 1| ≤ γ` entrywise with `n γ ≤ 1/2`, then `|‖c‖² - ‖Bᵀc‖²| ≤ 3 n γ ‖c‖²`. This is
   the plan's `Y (YᵀY)⁻¹ Yᵀ` lemma in the coordinates of `col Y`, where `Y (YᵀY)⁻¹ Yᵀ` is
   the identity and `Y Yᵀ` is `B Bᵀ`.
2. The eigenbasis expansion of the basis-free projector `specProj S T`
   (`specProj_apply_eq_sum`, `norm_sq_specProj_eq_sum`), its congruence in the eigenvalue
   set (`specProj_congr_of_iff`) and its additivity over disjoint sets
   (`specProj_union_of_disjoint`, `norm_sq_specProj_union`).
3. The residual bound (`exists_eigenvalue_near`, item (1) of the task) and its
   off-threshold form `sum_sq_off_le`: `∑_{λ_i ≤ τ} ⟨u_i, y⟩² ≤ (δ/m)²` for a vector `y`
   with `‖(S - ρ) y‖ ≤ δ` and `ρ ≥ τ + m`.
4. The two counts (plan 1.3 item 3 and item 4). `eigenvalues₀_le_of_split`: for
   `S = W + Q Qᵀ` with `Q` of `r` columns and `lamMax W ≤ τ`, every sorted eigenvalue of
   index `r` or more is at most `τ`. `le_card_filter_of_frame`: an `s`-frame of approximate
   eigenvectors above `τ + m` forces at least `s` eigenvalue indices above `τ`.
5. The frame lemma. `frame_core` is the coordinate form; `specProj_frame_approx` reads:
   for `s` approximate eigenvectors `y_k` (`‖(S - ρ_k) y_k‖ ≤ δ`, `ρ_k ≥ τ + m`, Gram
   entrywise within `δ'` of `1`, `s (δ' + δ/m) ≤ 1/2`) and exactly `s` eigenvalue indices
   above `τ`, `|‖P_{>τ} x‖² - ∑_k ⟨y_k, x⟩²| ≤ 5 s (δ' + δ/m) ‖x‖²` for every `x`.
   `specProjTop_frame_approx` is the same with `specProjTop S r` when the `r` indices
   above `τ` are the top `r`. Ties among the `ρ_k` are allowed; no block is named.
6. The split (plan 1.3 item 5). `specProjTop_eq_specProj_add_edge`: when every sorted
   eigenvalue of index `r` or more is at most `τ`, `specProjTop S r` is the sum of the
   projector on the eigenvalues above `τ` and the projector on the top-`r` eigenvalues at
   most `τ` (`P_edge`), and the squared norms add (`norm_sq_specProjTop_split`).

## Conventions

Vectors of the frame are `EuclideanSpace ℝ (Fin p)`; the coordinate lemmas use
`Fin p → ℝ` and `⬝ᵥ`. Gram hypotheses are entrywise, so every constant carries a factor
`s` or `s²` instead of an operator norm. The constant `5` is not sharp; the numeric check
`notes/archive/agent_reports/u3_frame.md` measured at most `0.12` of it (two seeds).
-/

open Matrix Finset
open scoped InnerProductSpace Matrix

namespace StackedSVD
namespace Frame

/-! ### 1. Entrywise bounds and the square frame lemma -/

section Elementary

/-- `(∑ |a_k|)² ≤ n ∑ a_k²`. -/
theorem sq_sum_abs_le {ι : Type*} [Fintype ι] (a : ι → ℝ) :
    (∑ k, |a k|) ^ 2 ≤ Fintype.card ι * ∑ k, a k ^ 2 := by
  have h := sq_sum_le_card_mul_sum_sq (s := Finset.univ) (f := fun k => |a k|)
  simpa [sq_abs, Finset.card_univ] using h

/-- An entrywise bound on `E` bounds the quadratic form `a ⬝ᵥ E a` by `n γ ‖a‖²`. -/
theorem abs_dotProduct_mulVec_le {n : ℕ} {E : Matrix (Fin n) (Fin n) ℝ} {γ : ℝ} (hγ : 0 ≤ γ)
    (hE : ∀ k l, |E k l| ≤ γ) (a : Fin n → ℝ) :
    |a ⬝ᵥ (E *ᵥ a)| ≤ n * γ * ∑ k, a k ^ 2 := by
  have h1 : |a ⬝ᵥ (E *ᵥ a)| ≤ ∑ k, ∑ l, γ * (|a k| * |a l|) := by
    simp only [dotProduct, Matrix.mulVec]
    calc |∑ k, a k * ∑ l, E k l * a l| ≤ ∑ k, |a k * ∑ l, E k l * a l| :=
          Finset.abs_sum_le_sum_abs _ _
      _ ≤ ∑ k, ∑ l, γ * (|a k| * |a l|) := by
          refine Finset.sum_le_sum fun k _ => ?_
          rw [abs_mul]
          calc |a k| * |∑ l, E k l * a l| ≤ |a k| * ∑ l, |E k l * a l| :=
                mul_le_mul_of_nonneg_left (Finset.abs_sum_le_sum_abs _ _) (abs_nonneg _)
            _ = ∑ l, |a k| * |E k l * a l| := Finset.mul_sum _ _ _
            _ ≤ ∑ l, γ * (|a k| * |a l|) := by
                refine Finset.sum_le_sum fun l _ => ?_
                rw [abs_mul]
                have := hE k l
                have := abs_nonneg (a k)
                have := abs_nonneg (a l)
                nlinarith [mul_nonneg (abs_nonneg (a k)) (abs_nonneg (a l))]
  have h2 : ∑ k, ∑ l, γ * (|a k| * |a l|) = γ * (∑ k, |a k|) ^ 2 := by
    rw [sq, Finset.sum_mul_sum]
    simp only [Finset.mul_sum]
  rw [h2] at h1
  refine h1.trans ?_
  have h3 := sq_sum_abs_le a
  rw [Fintype.card_fin] at h3
  calc γ * (∑ k, |a k|) ^ 2 ≤ γ * (n * ∑ k, a k ^ 2) := mul_le_mul_of_nonneg_left h3 hγ
    _ = n * γ * ∑ k, a k ^ 2 := by ring

/-- An entrywise bound on `E` bounds `‖E a‖² ≤ (n γ)² ‖a‖²`. -/
theorem sum_sq_mulVec_le {n : ℕ} {E : Matrix (Fin n) (Fin n) ℝ} {γ : ℝ}
    (hE : ∀ k l, |E k l| ≤ γ) (a : Fin n → ℝ) :
    ∑ k, (E *ᵥ a) k ^ 2 ≤ (n * γ) ^ 2 * ∑ k, a k ^ 2 := by
  have hk : ∀ k, |(E *ᵥ a) k| ≤ γ * ∑ l, |a l| := by
    intro k
    simp only [Matrix.mulVec, dotProduct]
    calc |∑ l, E k l * a l| ≤ ∑ l, |E k l * a l| := Finset.abs_sum_le_sum_abs _ _
      _ ≤ ∑ l, γ * |a l| := by
          refine Finset.sum_le_sum fun l _ => ?_
          rw [abs_mul]
          exact mul_le_mul_of_nonneg_right (hE k l) (abs_nonneg _)
      _ = γ * ∑ l, |a l| := (Finset.mul_sum _ _ _).symm
  have h3 := sq_sum_abs_le a
  rw [Fintype.card_fin] at h3
  have hk2 : ∀ k, (E *ᵥ a) k ^ 2 ≤ γ ^ 2 * (n * ∑ l, a l ^ 2) := by
    intro k
    rw [← sq_abs]
    calc |(E *ᵥ a) k| ^ 2 ≤ (γ * ∑ l, |a l|) ^ 2 :=
          pow_le_pow_left₀ (abs_nonneg _) (hk k) 2
      _ = γ ^ 2 * (∑ l, |a l|) ^ 2 := by ring
      _ ≤ γ ^ 2 * (n * ∑ l, a l ^ 2) := mul_le_mul_of_nonneg_left h3 (sq_nonneg _)
  calc ∑ k, (E *ᵥ a) k ^ 2 ≤ ∑ _k : Fin n, γ ^ 2 * (n * ∑ l, a l ^ 2) :=
        Finset.sum_le_sum fun k _ => hk2 k
    _ = (n * γ) ^ 2 * ∑ k, a k ^ 2 := by
        rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
        ring

/-- The Gram matrix `Bᵀ B` of the columns of `B`, entrywise. -/
theorem transpose_mul_apply {ι : Type*} [Fintype ι] {n : ℕ} (B : Matrix ι (Fin n) ℝ)
    (k l : Fin n) : (Bᵀ * B) k l = ∑ i, B i k * B i l := by
  simp [Matrix.mul_apply]

/-- `‖B a‖² = a ⬝ᵥ (BᵀB) a`. -/
theorem sum_sq_mulVec_eq {ι : Type*} [Fintype ι] {n : ℕ} (B : Matrix ι (Fin n) ℝ)
    (a : Fin n → ℝ) : ∑ i, (B *ᵥ a) i ^ 2 = a ⬝ᵥ ((Bᵀ * B) *ᵥ a) := by
  have h : (B *ᵥ a) ⬝ᵥ (B *ᵥ a) = a ⬝ᵥ ((Bᵀ * B) *ᵥ a) := by
    rw [← Matrix.mulVec_mulVec, Matrix.mulVec_transpose, Matrix.dotProduct_mulVec,
      dotProduct_comm]
  rw [← h]
  simp [dotProduct, sq]

/-- If `|BᵀB - 1| ≤ γ` entrywise with `n γ < 1`, then `B a = 0` forces `a = 0`. -/
theorem mulVec_eq_zero_of_gram_close {ι : Type*} [Fintype ι] {n : ℕ} {B : Matrix ι (Fin n) ℝ}
    {γ : ℝ} (hγ : 0 ≤ γ) (hB : ∀ k l, |(Bᵀ * B) k l - if k = l then 1 else 0| ≤ γ)
    (hsmall : n * γ < 1) {a : Fin n → ℝ} (ha : B *ᵥ a = 0) : a = 0 := by
  set E : Matrix (Fin n) (Fin n) ℝ := Bᵀ * B - 1 with hEdef
  have hE : ∀ k l, |E k l| ≤ γ := by
    intro k l
    simpa [hEdef, Matrix.one_apply] using hB k l
  have hBB : Bᵀ * B = 1 + E := by rw [hEdef]; abel
  have h0 : ∑ i, (B *ᵥ a) i ^ 2 = 0 := by simp [ha]
  rw [sum_sq_mulVec_eq, hBB, Matrix.add_mulVec, Matrix.one_mulVec, dotProduct_add] at h0
  have hsq : a ⬝ᵥ a = ∑ k, a k ^ 2 := by simp [dotProduct, sq]
  rw [hsq] at h0
  have h1 := abs_dotProduct_mulVec_le hγ hE a
  have hnn : 0 ≤ ∑ k, a k ^ 2 := Finset.sum_nonneg fun k _ => sq_nonneg _
  have hzero : ∑ k, a k ^ 2 = 0 := by
    rcases abs_le.mp h1 with ⟨h1l, _⟩
    nlinarith
  funext k
  have := (Finset.sum_eq_zero_iff_of_nonneg fun k _ => sq_nonneg (a k)).mp hzero k
    (Finset.mem_univ _)
  exact pow_eq_zero_iff (n := 2) (by norm_num) |>.mp this

/-- The columns of a matrix with `|BᵀB - 1| ≤ γ` entrywise and `n γ < 1` are linearly
independent. -/
theorem linearIndependent_of_gram_close {ι : Type*} [Fintype ι] {n : ℕ}
    {B : Matrix ι (Fin n) ℝ} {γ : ℝ} (hγ : 0 ≤ γ)
    (hB : ∀ k l, |(Bᵀ * B) k l - if k = l then 1 else 0| ≤ γ) (hsmall : n * γ < 1) :
    LinearIndependent ℝ (fun k : Fin n => fun i : ι => B i k) := by
  rw [Fintype.linearIndependent_iff]
  intro g hg
  have hBg : B *ᵥ g = 0 := by
    rw [← hg]
    funext i
    simp [Matrix.mulVec, dotProduct, mul_comm]
  have := mulVec_eq_zero_of_gram_close hγ hB hsmall hBg
  intro k
  rw [this]
  rfl

/-- A matrix with `n` nearly orthonormal columns has at least `n` rows. -/
theorem card_le_of_gram_close {ι : Type*} [Fintype ι] {n : ℕ}
    {B : Matrix ι (Fin n) ℝ} {γ : ℝ} (hγ : 0 ≤ γ)
    (hB : ∀ k l, |(Bᵀ * B) k l - if k = l then 1 else 0| ≤ γ) (hsmall : n * γ < 1) :
    n ≤ Fintype.card ι := by
  have h := (linearIndependent_of_gram_close hγ hB hsmall).fintype_card_le_finrank
  rwa [Fintype.card_fin, Module.finrank_fintype_fun_eq_card] at h

/-- The square frame lemma. For `B : Matrix ι (Fin n) ℝ` with `card ι = n` and
`|BᵀB - 1| ≤ γ` entrywise, `n γ ≤ 1/2`: `|‖c‖² - ‖Bᵀ c‖²| ≤ 3 n γ ‖c‖²` for every `c`.
In the plan's language (1.3 item 4): the projector `Y (YᵀY)⁻¹ Yᵀ` onto `col Y` is within
`3 n γ` of `Y Yᵀ` when `‖YᵀY - 1‖ ≤ γ`, read in the coordinates of `col Y`. -/
theorem sq_frame_close {ι : Type*} [Fintype ι] {n : ℕ} (hcard : Fintype.card ι = n)
    {B : Matrix ι (Fin n) ℝ} {γ : ℝ} (hγ : 0 ≤ γ)
    (hB : ∀ k l, |(Bᵀ * B) k l - if k = l then 1 else 0| ≤ γ) (hsmall : n * γ ≤ 1 / 2)
    (c : ι → ℝ) :
    |∑ i, c i ^ 2 - ∑ k, (Bᵀ *ᵥ c) k ^ 2| ≤ 3 * n * γ * ∑ i, c i ^ 2 := by
  set E : Matrix (Fin n) (Fin n) ℝ := Bᵀ * B - 1 with hEdef
  have hE : ∀ k l, |E k l| ≤ γ := by
    intro k l
    simpa [hEdef, Matrix.one_apply] using hB k l
  have hBB : Bᵀ * B = 1 + E := by rw [hEdef]; abel
  have hlt : n * γ < 1 := by linarith
  -- `B` is injective, hence surjective: `c = B a`
  have hinj : Function.Injective (Matrix.mulVecLin B) := by
    rw [← LinearMap.ker_eq_bot, LinearMap.ker_eq_bot']
    intro a ha
    exact mulVec_eq_zero_of_gram_close hγ hB hlt ha
  have hsurj : Function.Surjective (Matrix.mulVecLin B) := by
    have hfr : Module.finrank ℝ (Fin n → ℝ) = Module.finrank ℝ (ι → ℝ) := by
      rw [Module.finrank_fintype_fun_eq_card, Module.finrank_fintype_fun_eq_card, hcard,
        Fintype.card_fin]
    exact (LinearMap.injective_iff_surjective_of_finrank_eq_finrank hfr).mp hinj
  obtain ⟨a, ha⟩ := hsurj c
  rw [Matrix.mulVecLin_apply] at ha
  -- the three quantities
  have hBtc : Bᵀ *ᵥ c = a + E *ᵥ a := by
    rw [← ha, Matrix.mulVec_mulVec, hBB, Matrix.add_mulVec, Matrix.one_mulVec]
  have hC : ∑ i, c i ^ 2 = ∑ k, a k ^ 2 + a ⬝ᵥ (E *ᵥ a) := by
    rw [← ha, sum_sq_mulVec_eq, hBB, Matrix.add_mulVec, Matrix.one_mulVec, dotProduct_add]
    simp [dotProduct, sq]
  have hD : ∑ k, (Bᵀ *ᵥ c) k ^ 2
      = ∑ k, a k ^ 2 + 2 * a ⬝ᵥ (E *ᵥ a) + ∑ k, (E *ᵥ a) k ^ 2 := by
    rw [hBtc]
    simp only [Pi.add_apply, dotProduct]
    rw [Finset.mul_sum, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun k _ => ?_
    ring
  have h1 := abs_dotProduct_mulVec_le hγ hE a
  have h2 := sum_sq_mulVec_le hE a
  have hA : 0 ≤ ∑ k, a k ^ 2 := Finset.sum_nonneg fun k _ => sq_nonneg _
  have hu : 0 ≤ ∑ k, (E *ᵥ a) k ^ 2 := Finset.sum_nonneg fun k _ => sq_nonneg _
  have hnγ : 0 ≤ n * γ := by positivity
  rw [hC, hD]
  rcases abs_le.mp h1 with ⟨h1l, h1u⟩
  rw [abs_le]
  constructor <;> nlinarith [mul_le_mul_of_nonneg_left hsmall hnγ]

end Elementary


/-! ### 2. The eigenbasis expansion of `specProj` -/

section Expansion

variable {p : ℕ} {S : Matrix (Fin p) (Fin p) ℝ}

/-- The operator of a real symmetric matrix is symmetric. -/
theorem isSymmetric_toOp (hS : S.IsHermitian) : (toOp S).IsSymmetric :=
  Matrix.isSymmetric_toEuclideanLin_iff.mpr hS

/-- The eigenbasis diagonalizes the operator. -/
theorem toOp_eigvec (hS : S.IsHermitian) (i : Fin p) :
    toOp S (hS.eigenvectorBasis i) = hS.eigenvalues i • hS.eigenvectorBasis i := by
  apply WithLp.ofLp_injective
  simpa using hS.mulVec_eigenvectorBasis i

/-- The eigenbasis coordinate of `S x` is the eigenvalue times the coordinate of `x`. -/
theorem inner_eigvec_toOp (hS : S.IsHermitian) (x : EuclideanSpace ℝ (Fin p)) (i : Fin p) :
    ⟪hS.eigenvectorBasis i, toOp S x⟫_ℝ
      = hS.eigenvalues i * ⟪hS.eigenvectorBasis i, x⟫_ℝ := by
  rw [← isSymmetric_toOp hS, toOp_eigvec, real_inner_smul_left]

/-- An eigenvector of `S` for `t` is orthogonal to every basis eigenvector whose eigenvalue
differs from `t`. -/
theorem inner_eigvec_eq_zero_of_ne (hS : S.IsHermitian) {t : ℝ} {w : EuclideanSpace ℝ (Fin p)}
    (hw : w ∈ Module.End.eigenspace (toOp S) t) (i : Fin p) (hi : hS.eigenvalues i ≠ t) :
    ⟪hS.eigenvectorBasis i, w⟫_ℝ = 0 := by
  have hw' : toOp S w = t • w := Module.End.mem_eigenspace_iff.mp hw
  have h1 := inner_eigvec_toOp hS w i
  rw [hw', real_inner_smul_right] at h1
  have h2 : (hS.eigenvalues i - t) * ⟪hS.eigenvectorBasis i, w⟫_ℝ = 0 := by linarith
  rcases mul_eq_zero.mp h2 with h | h
  · exact absurd (sub_eq_zero.mp h) hi
  · exact h

/-- Parseval: `⟪x, y⟫ = ∑ ⟪u_i, x⟫ ⟪u_i, y⟫`. -/
theorem inner_eq_sum_coord (hS : S.IsHermitian) (x y : EuclideanSpace ℝ (Fin p)) :
    ⟪x, y⟫_ℝ = ∑ i, ⟪hS.eigenvectorBasis i, x⟫_ℝ * ⟪hS.eigenvectorBasis i, y⟫_ℝ := by
  rw [← hS.eigenvectorBasis.sum_inner_mul_inner x y]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [real_inner_comm x]

/-- Parseval: `‖x‖² = ∑ ⟪u_i, x⟫²`. -/
theorem norm_sq_eq_sum_coord (hS : S.IsHermitian) (x : EuclideanSpace ℝ (Fin p)) :
    ‖x‖ ^ 2 = ∑ i, ⟪hS.eigenvectorBasis i, x⟫_ℝ ^ 2 := by
  rw [← real_inner_self_eq_norm_sq, inner_eq_sum_coord hS]
  simp [sq]

/-- A basis eigenvector whose eigenvalue lies in `T` lies in `specSpace S T`. -/
theorem eigvec_mem_specSpace (hS : S.IsHermitian) {T : Set ℝ} (i : Fin p)
    (hi : hS.eigenvalues i ∈ T) : hS.eigenvectorBasis i ∈ specSpace S T :=
  Submodule.mem_iSup_of_mem (hS.eigenvalues i)
    (Submodule.mem_iSup_of_mem hi (Module.End.mem_eigenspace_iff.mpr (toOp_eigvec hS i)))

/-- `specProj S T` in the eigenbasis: keep the coordinates whose eigenvalue lies in `T`. -/
theorem specProj_apply_eq_sum (hS : S.IsHermitian) (T : Set ℝ) (x : EuclideanSpace ℝ (Fin p)) :
    specProj S T x = ∑ i, T.indicator (fun _ => ⟪hS.eigenvectorBasis i, x⟫_ℝ)
      (hS.eigenvalues i) • hS.eigenvectorBasis i := by
  classical
  set v : EuclideanSpace ℝ (Fin p) := ∑ i, T.indicator (fun _ => ⟪hS.eigenvectorBasis i, x⟫_ℝ)
      (hS.eigenvalues i) • hS.eigenvectorBasis i with hv
  have hproj : specProj S T x = (specSpace S T).starProjection x := rfl
  rw [hproj]
  refine Submodule.eq_starProjection_of_mem_of_inner_eq_zero ?_ ?_
  · refine Submodule.sum_mem _ fun i _ => ?_
    by_cases hi : hS.eigenvalues i ∈ T
    · exact Submodule.smul_mem _ _ (eigvec_mem_specSpace hS i hi)
    · rw [Set.indicator_of_notMem hi, zero_smul]
      exact Submodule.zero_mem _
  · intro w hw
    have hw' : w ∈ ⨆ (t : ℝ) (_ : t ∈ T), Module.End.eigenspace (toOp S) t := hw
    refine Submodule.iSup_induction (motive := fun w => ⟪x - v, w⟫_ℝ = 0) _ hw' ?_ (by simp) ?_
    · intro t w₁ hw₁
      refine Submodule.iSup_induction (motive := fun w => ⟪x - v, w⟫_ℝ = 0) _ hw₁ ?_ (by simp)
        ?_
      · intro ht w₂ hw₂
        rw [inner_sub_left, hv, sum_inner, inner_eq_sum_coord hS x w₂, ← Finset.sum_sub_distrib]
        refine Finset.sum_eq_zero fun i _ => ?_
        rw [real_inner_smul_left]
        by_cases hi : hS.eigenvalues i ∈ T
        · rw [Set.indicator_of_mem hi, sub_self]
        · rw [Set.indicator_of_notMem hi, zero_mul, sub_zero]
          have hne : hS.eigenvalues i ≠ t := fun h => hi (h ▸ ht)
          rw [inner_eigvec_eq_zero_of_ne hS hw₂ i hne, mul_zero]
      · intro u u' hu hu'
        rw [inner_add_right, hu, hu', add_zero]
    · intro u u' hu hu'
      rw [inner_add_right, hu, hu', add_zero]

/-- `‖specProj S T x‖² = ∑_{λ_i ∈ T} ⟪u_i, x⟫²`. -/
theorem norm_sq_specProj_eq_sum (hS : S.IsHermitian) (T : Set ℝ)
    (x : EuclideanSpace ℝ (Fin p)) :
    ‖specProj S T x‖ ^ 2
      = ∑ i, T.indicator (fun _ => ⟪hS.eigenvectorBasis i, x⟫_ℝ ^ 2) (hS.eigenvalues i) := by
  rw [specProj_apply_eq_sum hS, ← real_inner_self_eq_norm_sq,
    hS.eigenvectorBasis.orthonormal.inner_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  by_cases hi : hS.eigenvalues i ∈ T
  · simp [Set.indicator_of_mem hi, sq]
  · simp [Set.indicator_of_notMem hi]

/-- Two eigenvalue sets that agree on the eigenvalues give the same projector. -/
theorem specProj_congr_of_iff (hS : S.IsHermitian) {T T' : Set ℝ}
    (h : ∀ i, hS.eigenvalues i ∈ T ↔ hS.eigenvalues i ∈ T') :
    specProj S T = specProj S T' := by
  refine ContinuousLinearMap.ext fun x => ?_
  rw [specProj_apply_eq_sum hS, specProj_apply_eq_sum hS]
  refine Finset.sum_congr rfl fun i _ => ?_
  by_cases hi : hS.eigenvalues i ∈ T
  · rw [Set.indicator_of_mem hi, Set.indicator_of_mem ((h i).mp hi)]
  · rw [Set.indicator_of_notMem hi, Set.indicator_of_notMem fun h' => hi ((h i).mpr h')]

/-- The projector is additive over disjoint eigenvalue sets. -/
theorem specProj_union_of_disjoint (hS : S.IsHermitian) {T₁ T₂ : Set ℝ}
    (h : Disjoint T₁ T₂) :
    specProj S (T₁ ∪ T₂) = specProj S T₁ + specProj S T₂ := by
  refine ContinuousLinearMap.ext fun x => ?_
  rw [_root_.add_apply, specProj_apply_eq_sum hS, specProj_apply_eq_sum hS,
    specProj_apply_eq_sum hS, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [← add_smul, Set.indicator_union_of_disjoint h]

/-- Pythagoras for the projectors on two disjoint eigenvalue sets. -/
theorem norm_sq_specProj_union (hS : S.IsHermitian) {T₁ T₂ : Set ℝ} (h : Disjoint T₁ T₂)
    (x : EuclideanSpace ℝ (Fin p)) :
    ‖specProj S (T₁ ∪ T₂) x‖ ^ 2 = ‖specProj S T₁ x‖ ^ 2 + ‖specProj S T₂ x‖ ^ 2 := by
  rw [norm_sq_specProj_eq_sum hS, norm_sq_specProj_eq_sum hS, norm_sq_specProj_eq_sum hS,
    ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Set.indicator_union_of_disjoint h]

/-- The real inner product of `EuclideanSpace` is the dot product. -/
theorem inner_eq_dot (x y : EuclideanSpace ℝ (Fin p)) :
    ⟪x, y⟫_ℝ = WithLp.ofLp x ⬝ᵥ WithLp.ofLp y := by
  simp [EuclideanSpace.inner_eq_star_dotProduct, dotProduct_comm]

end Expansion

/-! ### 3. Residual bounds -/

section Residual

variable {p : ℕ} {S : Matrix (Fin p) (Fin p) ℝ}

/-- The residual in coordinates: `‖S y - ρ y‖² = ∑ (λ_i - ρ)² ⟪u_i, y⟫²`. -/
theorem norm_sq_residual (hS : S.IsHermitian) (ρ : ℝ) (y : EuclideanSpace ℝ (Fin p)) :
    ‖toOp S y - ρ • y‖ ^ 2
      = ∑ i, (hS.eigenvalues i - ρ) ^ 2 * ⟪hS.eigenvectorBasis i, y⟫_ℝ ^ 2 := by
  rw [norm_sq_eq_sum_coord hS]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [inner_sub_right, inner_eigvec_toOp, real_inner_smul_right]
  ring

/-- Item (1) of the task: a unit vector with residual `‖S y - ρ y‖ ≤ δ` puts an eigenvalue
of `S` within `δ` of `ρ`. -/
theorem exists_eigenvalue_near (hS : S.IsHermitian) {ρ δ : ℝ} {y : EuclideanSpace ℝ (Fin p)}
    (hy : ‖y‖ = 1) (hres : ‖toOp S y - ρ • y‖ ≤ δ) : ∃ i, |hS.eigenvalues i - ρ| ≤ δ := by
  by_contra h
  push Not at h
  have hδ : 0 ≤ δ := le_trans (norm_nonneg _) hres
  have hsq : ‖toOp S y - ρ • y‖ ^ 2 ≤ δ ^ 2 := pow_le_pow_left₀ (norm_nonneg _) hres 2
  rw [norm_sq_residual hS] at hsq
  have hone : ∑ i, ⟪hS.eigenvectorBasis i, y⟫_ℝ ^ 2 = 1 := by
    rw [← norm_sq_eq_sum_coord hS, hy, one_pow]
  obtain ⟨i₀, _, hi₀⟩ := Finset.exists_ne_zero_of_sum_ne_zero
    (show ∑ i, ⟪hS.eigenvectorBasis i, y⟫_ℝ ^ 2 ≠ 0 by rw [hone]; exact one_ne_zero)
  have hlt : δ ^ 2 * ∑ i, ⟪hS.eigenvectorBasis i, y⟫_ℝ ^ 2
      < ∑ i, (hS.eigenvalues i - ρ) ^ 2 * ⟪hS.eigenvectorBasis i, y⟫_ℝ ^ 2 := by
    rw [Finset.mul_sum]
    refine Finset.sum_lt_sum (fun i _ => ?_) ⟨i₀, Finset.mem_univ _, ?_⟩
    · have : δ ^ 2 ≤ (hS.eigenvalues i - ρ) ^ 2 := by
        rw [← sq_abs (hS.eigenvalues i - ρ)]
        exact pow_le_pow_left₀ hδ (h i).le 2
      exact mul_le_mul_of_nonneg_right this (sq_nonneg _)
    · have h1 : δ ^ 2 < (hS.eigenvalues i₀ - ρ) ^ 2 := by
        rw [← sq_abs (hS.eigenvalues i₀ - ρ)]
        exact pow_lt_pow_left₀ (h i₀) hδ (by norm_num)
      exact mul_lt_mul_of_pos_right h1 ((sq_nonneg _).lt_of_ne hi₀.symm)
  rw [hone, mul_one] at hlt
  linarith

/-- The mass of an approximate eigenvector on the eigenvalues at most `τ`: if
`∑ (λ_i - ρ)² b_i² ≤ δ²` and `ρ ≥ τ + m`, then `∑_{λ_i ≤ τ} b_i² ≤ (δ/m)²`. -/
theorem sum_sq_off_le {lam : Fin p → ℝ} {τ m ρ δ : ℝ} (hm : 0 < m) (hρ : τ + m ≤ ρ)
    {b : Fin p → ℝ} (hres : ∑ i, (lam i - ρ) ^ 2 * b i ^ 2 ≤ δ ^ 2) :
    ∑ i ∈ Finset.univ.filter (fun i => ¬ τ < lam i), b i ^ 2 ≤ (δ / m) ^ 2 := by
  have h1 : ∀ i ∈ Finset.univ.filter (fun i => ¬ τ < lam i),
      m ^ 2 * b i ^ 2 ≤ (lam i - ρ) ^ 2 * b i ^ 2 := by
    intro i hi
    rw [Finset.mem_filter] at hi
    have hle : lam i ≤ τ := not_lt.mp hi.2
    have hml : m ≤ ρ - lam i := by linarith
    have hm2 : m ^ 2 ≤ (lam i - ρ) ^ 2 := by
      nlinarith [mul_nonneg (sub_nonneg.mpr hml) (by linarith : 0 ≤ ρ - lam i + m)]
    exact mul_le_mul_of_nonneg_right hm2 (sq_nonneg _)
  have h2 : m ^ 2 * ∑ i ∈ Finset.univ.filter (fun i => ¬ τ < lam i), b i ^ 2
      ≤ ∑ i ∈ Finset.univ.filter (fun i => ¬ τ < lam i), (lam i - ρ) ^ 2 * b i ^ 2 := by
    rw [Finset.mul_sum]
    exact Finset.sum_le_sum h1
  have h3 : ∑ i ∈ Finset.univ.filter (fun i => ¬ τ < lam i), (lam i - ρ) ^ 2 * b i ^ 2
      ≤ ∑ i, (lam i - ρ) ^ 2 * b i ^ 2 :=
    Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _)
      fun i _ _ => by positivity
  rw [div_pow, le_div_iff₀ (by positivity)]
  linarith

end Residual

/-! ### 4. Counts: at most `r` indices above `lamMax W`, at least `s` above a frame -/

section Count

variable {p : ℕ} {S : Matrix (Fin p) (Fin p) ℝ}

/-- `0 ≤ v ⬝ᵥ v` for real vectors. -/
theorem dot_self_nonneg {n : ℕ} (v : Fin n → ℝ) : 0 ≤ v ⬝ᵥ v :=
  Finset.sum_nonneg fun i _ => mul_self_nonneg (v i)

/-- Plan 1.3 item 3. For `S = W + Q Qᵀ` with `Q` of `r` columns and `lamMax W ≤ τ`, every
sorted eigenvalue of `S` with index at least `r` is at most `τ`. The eigenvectors of the
top `k + 1` eigenvalues span a space on which `Qᵀ` is injective, so `k + 1 ≤ r`. -/
theorem eigenvalues₀_le_of_split {W : Matrix (Fin p) (Fin p) ℝ} {r : ℕ}
    {Q : Matrix (Fin p) (Fin r) ℝ} (hW : W.IsHermitian) (hS : S.IsHermitian)
    (hSeq : S = W + Q * Qᵀ) {τ : ℝ} (hτ : lamMax W hW ≤ τ)
    (k : Fin (Fintype.card (Fin p))) (hk : r ≤ (k : ℕ)) : hS.eigenvalues₀ k ≤ τ := by
  by_contra hlt
  push Not at hlt
  have hnp : (k : ℕ) + 1 ≤ Fintype.card (Fin p) := k.isLt
  set u : Fin ((k : ℕ) + 1) → EuclideanSpace ℝ (Fin p) :=
    fun j => hS.eigenvectorBasis (eigIdx p (Fin.castLE hnp j)) with hu
  have hortho : Orthonormal ℝ u :=
    hS.eigenvectorBasis.orthonormal.comp _
      ((eigIdx p).injective.comp (Fin.castLE_injective hnp))
  set lam : Fin ((k : ℕ) + 1) → ℝ := fun j => hS.eigenvalues (eigIdx p (Fin.castLE hnp j))
    with hlam
  have heig : ∀ j, τ < lam j := by
    intro j
    rw [hlam]
    simp only
    rw [eigenvalues_eigIdx]
    refine lt_of_lt_of_le hlt (hS.eigenvalues₀_antitone ?_)
    rw [Fin.le_def, Fin.val_castLE]
    omega
  have hli : LinearIndependent ℝ fun j => Qᵀ *ᵥ WithLp.ofLp (u j) := by
    rw [Fintype.linearIndependent_iff]
    intro g hg
    set x : EuclideanSpace ℝ (Fin p) := ∑ j, g j • u j with hx
    have hQx : Qᵀ *ᵥ WithLp.ofLp x = 0 := by
      rw [← hg, hx, WithLp.ofLp_sum, Matrix.mulVec_sum]
      refine Finset.sum_congr rfl fun j _ => ?_
      rw [WithLp.ofLp_smul, Matrix.mulVec_smul]
    have hSx : toOp S x = ∑ j, (g j * lam j) • u j := by
      rw [hx, map_sum]
      refine Finset.sum_congr rfl fun j _ => ?_
      rw [map_smul, hu, hlam]
      simp only
      rw [toOp_eigvec, smul_smul]
    have hq1 : ⟪x, toOp S x⟫_ℝ = ∑ j, g j * (g j * lam j) := by
      rw [hSx, hx, hortho.inner_sum]
      simp
    have hq0 : ⟪x, x⟫_ℝ = ∑ j, g j * g j := by
      rw [hx, hortho.inner_sum]
      simp
    have hq2 : ⟪x, toOp S x⟫_ℝ ≤ τ * ⟪x, x⟫_ℝ := by
      have hofLp : WithLp.ofLp (toOp S x) = S *ᵥ WithLp.ofLp x := rfl
      rw [inner_eq_dot, inner_eq_dot, hofLp, hSeq, Matrix.add_mulVec, dotProduct_add,
        ← Matrix.mulVec_mulVec, hQx, Matrix.mulVec_zero, dotProduct_zero, add_zero]
      calc WithLp.ofLp x ⬝ᵥ (W *ᵥ WithLp.ofLp x)
          ≤ lamMax W hW * (WithLp.ofLp x ⬝ᵥ WithLp.ofLp x) :=
            R4.dotProduct_mulVec_le_lamMax hW _
        _ ≤ τ * (WithLp.ofLp x ⬝ᵥ WithLp.ofLp x) :=
            mul_le_mul_of_nonneg_right hτ (dot_self_nonneg _)
    have hsum : ∑ j, g j * g j * (lam j - τ) ≤ 0 := by
      have heq : ∑ j, g j * g j * (lam j - τ)
          = ∑ j, g j * (g j * lam j) - τ * ∑ j, g j * g j := by
        rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
        refine Finset.sum_congr rfl fun j _ => ?_
        ring
      rw [heq]
      rw [hq1, hq0] at hq2
      linarith
    intro j
    have hnn : ∀ j' ∈ Finset.univ, 0 ≤ g j' * g j' * (lam j' - τ) := fun j' _ =>
      mul_nonneg (mul_self_nonneg _) (sub_nonneg.mpr (heig j').le)
    have hj : g j * g j * (lam j - τ) ≤ 0 :=
      le_trans (Finset.single_le_sum hnn (Finset.mem_univ j)) hsum
    have hpos : 0 < lam j - τ := sub_pos.mpr (heig j)
    have hgg : g j * g j ≤ 0 := by nlinarith [mul_self_nonneg (g j)]
    exact mul_self_eq_zero.mp (le_antisymm hgg (mul_self_nonneg _))
  have hcard := hli.fintype_card_le_finrank
  rw [Module.finrank_fintype_fun_eq_card] at hcard
  simp only [Fintype.card_fin] at hcard
  omega

/-- The count of eigenvalue indices with a property, in the matrix index and in the sorted
index. -/
theorem card_filter_eigenvalues (hS : S.IsHermitian) (P : ℝ → Prop) [DecidablePred P] :
    (Finset.univ.filter fun i => P (hS.eigenvalues i)).card
      = (Finset.univ.filter fun k => P (hS.eigenvalues₀ k)).card := by
  symm
  refine Finset.card_bij (fun k _ => eigIdx p k) ?_ ?_ ?_
  · intro k hk
    rw [Finset.mem_filter] at hk ⊢
    refine ⟨Finset.mem_univ _, ?_⟩
    rw [eigenvalues_eigIdx]
    exact hk.2
  · intro k₁ _ k₂ _ h
    exact (eigIdx p).injective h
  · intro i hi
    refine ⟨(eigIdx p).symm i, ?_, Equiv.apply_symm_apply _ _⟩
    rw [Finset.mem_filter] at hi ⊢
    refine ⟨Finset.mem_univ _, ?_⟩
    rw [← eigenvalues_eigIdx hS, Equiv.apply_symm_apply]
    exact hi.2

/-- The card form of `eigenvalues₀_le_of_split`: at most `r` eigenvalue indices of
`W + Q Qᵀ` lie above `lamMax W`. -/
theorem card_filter_le_of_split {W : Matrix (Fin p) (Fin p) ℝ} {r : ℕ}
    {Q : Matrix (Fin p) (Fin r) ℝ} (hW : W.IsHermitian) (hS : S.IsHermitian)
    (hSeq : S = W + Q * Qᵀ) {τ : ℝ} (hτ : lamMax W hW ≤ τ) :
    (Finset.univ.filter fun i => τ < hS.eigenvalues i).card ≤ r := by
  rw [card_filter_eigenvalues hS (fun t => τ < t)]
  calc (Finset.univ.filter fun k => τ < hS.eigenvalues₀ k).card
      ≤ (Finset.univ.filter fun k : Fin (Fintype.card (Fin p)) => (k : ℕ) < r).card := by
        refine Finset.card_le_card fun k hk => ?_
        rw [Finset.mem_filter] at hk ⊢
        refine ⟨hk.1, ?_⟩
        by_contra hge
        exact absurd hk.2
          (not_lt.mpr (eigenvalues₀_le_of_split hW hS hSeq hτ k (not_lt.mp hge)))
    _ ≤ (Finset.range r).card := by
        refine Finset.card_le_card_of_injOn (fun k => (k : ℕ)) ?_ ?_
        · intro k hk
          rw [Finset.mem_coe, Finset.mem_filter] at hk
          rw [Finset.mem_coe, Finset.mem_range]
          exact hk.2
        · intro a _ b _ h
          exact Fin.ext h
    _ = r := Finset.card_range r

/-- Exactly `r` sorted indices above `τ` gives the top-`r` gap. -/
theorem topGap_of_count (hS : S.IsHermitian) {r : ℕ} {τ : ℝ}
    (hcount : ∀ k, τ < hS.eigenvalues₀ k ↔ (k : ℕ) < r) : TopGap S hS r := by
  intro k l hk hl
  have h1 := (hcount k).mpr hk
  have h2 : ¬ τ < hS.eigenvalues₀ l := fun h => by
    have := (hcount l).mp h
    omega
  push Not at h2
  linarith

/-- Exactly `r` sorted indices above `τ` gives `r` matrix indices above `τ`. -/
theorem card_filter_of_count (hS : S.IsHermitian) {r : ℕ} (hrp : r ≤ p) {τ : ℝ}
    (hcount : ∀ k, τ < hS.eigenvalues₀ k ↔ (k : ℕ) < r) :
    (Finset.univ.filter fun i => τ < hS.eigenvalues i).card = r := by
  rw [card_filter_eigenvalues hS (fun t => τ < t)]
  have hset : (Finset.univ.filter fun k => τ < hS.eigenvalues₀ k)
      = Finset.univ.filter fun k : Fin (Fintype.card (Fin p)) => (k : ℕ) < r := by
    ext k
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    exact hcount k
  rw [hset, card_filter_lt (by simpa using hrp)]

/-- When the `r` sorted indices above `τ` are exactly the top `r`, `specProjTop S r` is the
projector on the eigenvalues above `τ`. -/
theorem specProjTop_eq_specProj_Ioi (hS : S.IsHermitian) {r : ℕ} {τ : ℝ}
    (hcount : ∀ k, τ < hS.eigenvalues₀ k ↔ (k : ℕ) < r) :
    specProjTop S hS r = specProj S (Set.Ioi τ) := by
  have hgap := topGap_of_count hS hcount
  change specProj S (topEigSet S hS r) = specProj S (Set.Ioi τ)
  refine specProj_congr_of_iff hS fun i => ?_
  have hi : hS.eigenvalues i = hS.eigenvalues₀ ((eigIdx p).symm i) := by
    rw [← eigenvalues_eigIdx hS, Equiv.apply_symm_apply]
  rw [hi, mem_topEigSet_iff hS hgap, Set.mem_Ioi, hcount]

end Count

/-! ### 5. The frame lemma -/

section FrameLemma

variable {p : ℕ}

/-- The Gram entry of the restrictions to `I = {i | τ < λ_i}` is within `δ' + (δ/m)²` of
`δ_kl`. -/
theorem abs_gram_on_sub_le {lam : Fin p → ℝ} {τ m δ δ' : ℝ} {s : ℕ} (hm : 0 < m)
    {b : Fin s → Fin p → ℝ} {ρ : Fin s → ℝ} (hρ : ∀ k, τ + m ≤ ρ k)
    (hres : ∀ k, ∑ i, (lam i - ρ k) ^ 2 * b k i ^ 2 ≤ δ ^ 2)
    (hgram : ∀ k l, |b k ⬝ᵥ b l - if k = l then 1 else 0| ≤ δ') (k l : Fin s) :
    |∑ i ∈ Finset.univ.filter (fun i => τ < lam i), b k i * b l i - if k = l then 1 else 0|
      ≤ δ' + (δ / m) ^ 2 := by
  have hsplit := Finset.sum_filter_add_sum_filter_not Finset.univ (fun i => τ < lam i)
    (fun i => b k i * b l i)
  have hoff : |∑ i ∈ Finset.univ.filter (fun i => ¬ τ < lam i), b k i * b l i|
      ≤ (δ / m) ^ 2 := by
    have hk := sum_sq_off_le hm (hρ k) (hres k)
    have hl := sum_sq_off_le hm (hρ l) (hres l)
    refine abs_le_of_sq_le_sq ?_ (by positivity)
    calc (∑ i ∈ Finset.univ.filter (fun i => ¬ τ < lam i), b k i * b l i) ^ 2
        ≤ (∑ i ∈ Finset.univ.filter (fun i => ¬ τ < lam i), b k i ^ 2)
          * ∑ i ∈ Finset.univ.filter (fun i => ¬ τ < lam i), b l i ^ 2 :=
          Finset.sum_mul_sq_le_sq_mul_sq _ _ _
      _ ≤ (δ / m) ^ 2 * (δ / m) ^ 2 :=
          mul_le_mul hk hl (Finset.sum_nonneg fun i _ => sq_nonneg _) (by positivity)
      _ = ((δ / m) ^ 2) ^ 2 := by ring
  have hdot : b k ⬝ᵥ b l = ∑ i, b k i * b l i := rfl
  have hg := hgram k l
  rw [hdot, ← hsplit] at hg
  obtain ⟨hg1, hg2⟩ := abs_le.mp hg
  obtain ⟨ho1, ho2⟩ := abs_le.mp hoff
  rw [abs_le]
  constructor <;> linarith

/-- Plan 1.3 item 4, the rank lower bound: an `s`-frame of approximate eigenvectors above
`τ + m` forces at least `s` eigenvalue indices above `τ`. -/
theorem le_card_filter_of_frame {lam : Fin p → ℝ} {τ m δ δ' : ℝ} {s : ℕ} (hm : 0 < m)
    (hδ' : 0 ≤ δ') {b : Fin s → Fin p → ℝ} {ρ : Fin s → ℝ} (hρ : ∀ k, τ + m ≤ ρ k)
    (hres : ∀ k, ∑ i, (lam i - ρ k) ^ 2 * b k i ^ 2 ≤ δ ^ 2)
    (hgram : ∀ k l, |b k ⬝ᵥ b l - if k = l then 1 else 0| ≤ δ')
    (hsmall : s * (δ' + (δ / m) ^ 2) < 1) :
    s ≤ (Finset.univ.filter fun i => τ < lam i).card := by
  set I := Finset.univ.filter fun i => τ < lam i with hI
  let B : Matrix I (Fin s) ℝ := fun i k => b k i
  have hBB : ∀ k l, |(Bᵀ * B) k l - if k = l then 1 else 0| ≤ δ' + (δ / m) ^ 2 := by
    intro k l
    rw [transpose_mul_apply]
    have hsum : ∑ i : I, B i k * B i l = ∑ i ∈ I, b k i * b l i :=
      Finset.sum_coe_sort I (fun i => b k i * b l i)
    rw [hsum]
    exact abs_gram_on_sub_le hm hρ hres hgram k l
  have h := card_le_of_gram_close (by positivity) hBB hsmall
  rwa [Fintype.card_coe] at h

/-- The cross term of the frame lemma: `|(α + β)² - α²| ≤ (η (2 + δ') + η²) C` when
`α² ≤ (1 + δ') C` and `β² ≤ η² C`. -/
theorem abs_sq_add_sub_sq_le {α β η δ' C : ℝ} (hη : 0 ≤ η)
    (hα : α ^ 2 ≤ (1 + δ') * C) (hβ : β ^ 2 ≤ η ^ 2 * C) :
    |(α + β) ^ 2 - α ^ 2| ≤ (η * (2 + δ') + η ^ 2) * C := by
  have hβ0 : 0 ≤ β ^ 2 := sq_nonneg _
  have habs : α * β ≤ |α| * |β| := by rw [← abs_mul]; exact le_abs_self _
  have habs' : -(α * β) ≤ |α| * |β| := by rw [← abs_mul]; exact neg_le_abs _
  have hcross : 2 * (|α| * |β|) ≤ η * (1 + δ') * C + η * C := by
    rcases eq_or_lt_of_le hη with hη0 | hηpos
    · subst hη0
      have hβz : β = 0 := by
        have : β ^ 2 ≤ 0 := by simpa using hβ
        exact pow_eq_zero_iff (n := 2) (by norm_num) |>.mp (le_antisymm this hβ0)
      simp [hβz]
    · have hkey : η * (2 * (|α| * |β|)) ≤ η * (η * (1 + δ') * C + η * C) := by
        have h1 : η ^ 2 * α ^ 2 ≤ η ^ 2 * ((1 + δ') * C) :=
          mul_le_mul_of_nonneg_left hα (sq_nonneg _)
        have hαa : |α| ^ 2 = α ^ 2 := sq_abs α
        have hβa : |β| ^ 2 = β ^ 2 := sq_abs β
        nlinarith [sq_nonneg (η * |α| - |β|)]
      exact le_of_mul_le_mul_left hkey hηpos
  have hexp : (α + β) ^ 2 - α ^ 2 = 2 * (α * β) + β ^ 2 := by ring
  rw [hexp, abs_le]
  constructor <;> nlinarith

/-- The frame lemma in the eigenbasis coordinates (plan 1.3 item 4). `lam` is the eigenvalue
list, `I = {i | τ < λ_i}` has `s` elements, `b_k` is the coordinate vector of the `k`-th
approximate eigenvector, and `c` that of the test vector. -/
theorem frame_core {lam : Fin p → ℝ} {τ m δ δ' : ℝ} {s : ℕ} (hm : 0 < m) (hδ : 0 ≤ δ)
    (hδ' : 0 ≤ δ') (hI : (Finset.univ.filter fun i => τ < lam i).card = s)
    {b : Fin s → Fin p → ℝ} {ρ : Fin s → ℝ} (hρ : ∀ k, τ + m ≤ ρ k)
    (hres : ∀ k, ∑ i, (lam i - ρ k) ^ 2 * b k i ^ 2 ≤ δ ^ 2)
    (hgram : ∀ k l, |b k ⬝ᵥ b l - if k = l then 1 else 0| ≤ δ')
    (hsmall : s * (δ' + δ / m) ≤ 1 / 2) (c : Fin p → ℝ) :
    |∑ i ∈ Finset.univ.filter (fun i => τ < lam i), c i ^ 2 - ∑ k, (b k ⬝ᵥ c) ^ 2|
      ≤ 5 * s * (δ' + δ / m) * ∑ i, c i ^ 2 := by
  set I := Finset.univ.filter fun i => τ < lam i with hI_def
  set η := δ / m with hη
  have hη0 : 0 ≤ η := by positivity
  rcases Nat.eq_zero_or_pos s with hs0 | hs
  · subst hs0
    have hIe : I = ∅ := Finset.card_eq_zero.mp hI
    simp [hIe]
  have hs1 : (1 : ℝ) ≤ s := by exact_mod_cast hs
  have hs0 : (0 : ℝ) ≤ s := by linarith
  have hsum_half : δ' + η ≤ 1 / 2 := by nlinarith
  have hδ'half : δ' ≤ 1 / 2 := by linarith
  have hηhalf : η ≤ 1 / 2 := by linarith
  set C := ∑ i, c i ^ 2 with hC
  have hC0 : 0 ≤ C := Finset.sum_nonneg fun i _ => sq_nonneg _
  set α : Fin s → ℝ := fun k => ∑ i ∈ I, b k i * c i with hα
  set β : Fin s → ℝ := fun k => ∑ i ∈ Finset.univ.filter (fun i => ¬ τ < lam i), b k i * c i
    with hβ
  have hsplit : ∀ k, b k ⬝ᵥ c = α k + β k := fun k =>
    (Finset.sum_filter_add_sum_filter_not Finset.univ (fun i => τ < lam i)
      (fun i => b k i * c i)).symm
  have hCI : ∑ i ∈ I, c i ^ 2 ≤ C :=
    Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _) fun i _ _ => sq_nonneg _
  have hCoff : ∑ i ∈ Finset.univ.filter (fun i => ¬ τ < lam i), c i ^ 2 ≤ C :=
    Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _) fun i _ _ => sq_nonneg _
  have hbk : ∀ k, ∑ i ∈ I, b k i ^ 2 ≤ 1 + δ' := by
    intro k
    have h1 : ∑ i ∈ I, b k i ^ 2 ≤ ∑ i, b k i ^ 2 :=
      Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _) fun i _ _ => sq_nonneg _
    have h2 := hgram k k
    have h3 : b k ⬝ᵥ b k = ∑ i, b k i ^ 2 := by simp [dotProduct, sq]
    rw [if_pos rfl, h3] at h2
    linarith [(abs_le.mp h2).2]
  have hα2 : ∀ k, α k ^ 2 ≤ (1 + δ') * C := by
    intro k
    calc α k ^ 2 ≤ (∑ i ∈ I, b k i ^ 2) * ∑ i ∈ I, c i ^ 2 :=
          Finset.sum_mul_sq_le_sq_mul_sq _ _ _
      _ ≤ (1 + δ') * C :=
          mul_le_mul (hbk k) hCI (Finset.sum_nonneg fun i _ => sq_nonneg _) (by linarith)
  have hβ2 : ∀ k, β k ^ 2 ≤ η ^ 2 * C := by
    intro k
    calc β k ^ 2 ≤ (∑ i ∈ Finset.univ.filter (fun i => ¬ τ < lam i), b k i ^ 2)
          * ∑ i ∈ Finset.univ.filter (fun i => ¬ τ < lam i), c i ^ 2 :=
          Finset.sum_mul_sq_le_sq_mul_sq _ _ _
      _ ≤ η ^ 2 * C :=
          mul_le_mul (sum_sq_off_le hm (hρ k) (hres k)) hCoff
            (Finset.sum_nonneg fun i _ => sq_nonneg _) (by positivity)
  -- the cross term, summed over `k`
  have hcross : |∑ k, (b k ⬝ᵥ c) ^ 2 - ∑ k, α k ^ 2| ≤ s * ((η * (2 + δ') + η ^ 2) * C) := by
    rw [← Finset.sum_sub_distrib]
    refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
    calc ∑ k, |(b k ⬝ᵥ c) ^ 2 - α k ^ 2|
        ≤ ∑ _k : Fin s, (η * (2 + δ') + η ^ 2) * C := by
          refine Finset.sum_le_sum fun k _ => ?_
          rw [hsplit k]
          exact abs_sq_add_sub_sq_le hη0 (hα2 k) (hβ2 k)
      _ = s * ((η * (2 + δ') + η ^ 2) * C) := by
          rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  -- the `I` part, by the square frame lemma
  set γ := δ' + η ^ 2 with hγ
  have hγ0 : 0 ≤ γ := by positivity
  have hη2 : η ^ 2 ≤ η := by nlinarith
  have hsγ : s * γ ≤ 1 / 2 := by nlinarith
  have hmain : |∑ i ∈ I, c i ^ 2 - ∑ k, α k ^ 2| ≤ 3 * s * γ * ∑ i ∈ I, c i ^ 2 := by
    let B : Matrix I (Fin s) ℝ := fun i k => b k i
    have hBB : ∀ k l, |(Bᵀ * B) k l - if k = l then 1 else 0| ≤ γ := by
      intro k l
      rw [transpose_mul_apply]
      have hsum : ∑ i : I, B i k * B i l = ∑ i ∈ I, b k i * b l i :=
        Finset.sum_coe_sort I (fun i => b k i * b l i)
      rw [hsum]
      exact abs_gram_on_sub_le hm hρ hres hgram k l
    have hcard : Fintype.card I = s := by rw [Fintype.card_coe, hI]
    have h := sq_frame_close hcard hγ0 hBB hsγ (fun i : I => c i)
    have hc1 : ∑ i : I, c i ^ 2 = ∑ i ∈ I, c i ^ 2 := Finset.sum_coe_sort I (fun i => c i ^ 2)
    have hc2 : ∀ k, (Bᵀ *ᵥ fun i : I => c i) k = α k := by
      intro k
      simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply]
      exact Finset.sum_coe_sort I (fun i => b k i * c i)
    simp only [hc2, hc1] at h
    exact h
  -- assemble
  have hconst : 3 * s * γ + s * (η * (2 + δ') + η ^ 2) ≤ 5 * s * (δ' + η) := by
    have hηδ : η * δ' ≤ δ' / 2 := by nlinarith
    nlinarith
  have hI_nonneg : 0 ≤ ∑ i ∈ I, c i ^ 2 := Finset.sum_nonneg fun i _ => sq_nonneg _
  have h3 : 3 * s * γ * ∑ i ∈ I, c i ^ 2 ≤ 3 * s * γ * C :=
    mul_le_mul_of_nonneg_left hCI (by positivity)
  calc |∑ i ∈ I, c i ^ 2 - ∑ k, (b k ⬝ᵥ c) ^ 2|
      ≤ |∑ i ∈ I, c i ^ 2 - ∑ k, α k ^ 2| + |∑ k, (b k ⬝ᵥ c) ^ 2 - ∑ k, α k ^ 2| := by
        rw [abs_le]
        constructor <;> linarith [abs_le.mp (le_refl |∑ i ∈ I, c i ^ 2 - ∑ k, α k ^ 2|),
          abs_le.mp (le_refl |∑ k, (b k ⬝ᵥ c) ^ 2 - ∑ k, α k ^ 2|)]
    _ ≤ 3 * s * γ * C + s * ((η * (2 + δ') + η ^ 2) * C) :=
        add_le_add (hmain.trans h3) hcross
    _ = (3 * s * γ + s * (η * (2 + δ') + η ^ 2)) * C := by ring
    _ ≤ 5 * s * (δ' + η) * C := mul_le_mul_of_nonneg_right hconst hC0

variable {S : Matrix (Fin p) (Fin p) ℝ}

/-- The frame lemma (plan 1.3 item 4, task item (4)). Let `S` have exactly `s` eigenvalue
indices above `τ`, and let `y_1, ..., y_s` satisfy `‖S y_k - ρ_k y_k‖ ≤ δ` with
`ρ_k ≥ τ + m`, Gram entries within `δ'` of `δ_kl`, and `s (δ' + δ/m) ≤ 1/2`. Then the
projector `P_{>τ}` on the eigenvalues above `τ` satisfies
`|‖P_{>τ} x‖² - ∑_k ⟪y_k, x⟫²| ≤ 5 s (δ' + δ/m) ‖x‖²` for every `x`. This is the quadratic
form of `‖P_{>τ} - ∑_k y_k y_kᵀ‖ ≤ 5 s (δ' + δ/m)`. Ties among the `ρ_k` are allowed. -/
theorem specProj_frame_approx (hS : S.IsHermitian) {τ m δ δ' : ℝ} {s : ℕ}
    (hm : 0 < m) (hδ : 0 ≤ δ) (hδ' : 0 ≤ δ')
    (hI : (Finset.univ.filter fun i => τ < hS.eigenvalues i).card = s)
    {y : Fin s → EuclideanSpace ℝ (Fin p)} {ρ : Fin s → ℝ} (hρ : ∀ k, τ + m ≤ ρ k)
    (hres : ∀ k, ‖toOp S (y k) - ρ k • y k‖ ≤ δ)
    (hgram : ∀ k l, |⟪y k, y l⟫_ℝ - if k = l then 1 else 0| ≤ δ')
    (hsmall : s * (δ' + δ / m) ≤ 1 / 2) (x : EuclideanSpace ℝ (Fin p)) :
    |‖specProj S (Set.Ioi τ) x‖ ^ 2 - ∑ k, ⟪y k, x⟫_ℝ ^ 2|
      ≤ 5 * s * (δ' + δ / m) * ‖x‖ ^ 2 := by
  set b : Fin s → Fin p → ℝ := fun k i => ⟪hS.eigenvectorBasis i, y k⟫_ℝ with hb
  set c : Fin p → ℝ := fun i => ⟪hS.eigenvectorBasis i, x⟫_ℝ with hc
  have hres' : ∀ k, ∑ i, (hS.eigenvalues i - ρ k) ^ 2 * b k i ^ 2 ≤ δ ^ 2 := by
    intro k
    have h := pow_le_pow_left₀ (norm_nonneg _) (hres k) 2
    rw [norm_sq_residual hS] at h
    exact h
  have hgram' : ∀ k l, |b k ⬝ᵥ b l - if k = l then 1 else 0| ≤ δ' := by
    intro k l
    have h : b k ⬝ᵥ b l = ⟪y k, y l⟫_ℝ := by
      rw [inner_eq_sum_coord hS]
      rfl
    rw [h]
    exact hgram k l
  have hP : ‖specProj S (Set.Ioi τ) x‖ ^ 2
      = ∑ i ∈ Finset.univ.filter (fun i => τ < hS.eigenvalues i), c i ^ 2 := by
    rw [norm_sq_specProj_eq_sum hS, Finset.sum_filter]
    refine Finset.sum_congr rfl fun i _ => ?_
    by_cases hi : τ < hS.eigenvalues i
    · rw [if_pos hi, Set.indicator_of_mem (Set.mem_Ioi.mpr hi)]
    · rw [if_neg hi, Set.indicator_of_notMem fun h => hi (Set.mem_Ioi.mp h)]
  have hyx : ∀ k, ⟪y k, x⟫_ℝ = b k ⬝ᵥ c := by
    intro k
    rw [inner_eq_sum_coord hS]
    rfl
  have hx : ‖x‖ ^ 2 = ∑ i, c i ^ 2 := norm_sq_eq_sum_coord hS x
  rw [hP, hx]
  simp only [hyx]
  exact frame_core hm hδ hδ' hI hρ hres' hgram' hsmall c

/-- The frame lemma for `specProjTop S r`: when the `r` sorted indices above `τ` are exactly
the top `r` (`hcount`), `|‖specProjTop S r x‖² - ∑_k ⟪y_k, x⟫²| ≤ 5 r (δ' + δ/m) ‖x‖²`. -/
theorem specProjTop_frame_approx (hS : S.IsHermitian) {r : ℕ} (hrp : r ≤ p) {τ : ℝ}
    (hcount : ∀ k, τ < hS.eigenvalues₀ k ↔ (k : ℕ) < r) {m δ δ' : ℝ}
    (hm : 0 < m) (hδ : 0 ≤ δ) (hδ' : 0 ≤ δ')
    {y : Fin r → EuclideanSpace ℝ (Fin p)} {ρ : Fin r → ℝ} (hρ : ∀ k, τ + m ≤ ρ k)
    (hres : ∀ k, ‖toOp S (y k) - ρ k • y k‖ ≤ δ)
    (hgram : ∀ k l, |⟪y k, y l⟫_ℝ - if k = l then 1 else 0| ≤ δ')
    (hsmall : r * (δ' + δ / m) ≤ 1 / 2) (x : EuclideanSpace ℝ (Fin p)) :
    |‖specProjTop S hS r x‖ ^ 2 - ∑ k, ⟪y k, x⟫_ℝ ^ 2|
      ≤ 5 * r * (δ' + δ / m) * ‖x‖ ^ 2 := by
  rw [specProjTop_eq_specProj_Ioi hS hcount]
  exact specProj_frame_approx hS hm hδ hδ' (card_filter_of_count hS hrp hcount) hρ hres hgram
    hsmall x

/-- The frame lemma with `Fin p → ℝ` vectors and dot products, the form the resolvent
layer (`RMT/R4.lean`) produces. -/
theorem specProj_frame_approx_dot (hS : S.IsHermitian) {τ m δ δ' : ℝ} {s : ℕ}
    (hm : 0 < m) (hδ : 0 ≤ δ) (hδ' : 0 ≤ δ')
    (hI : (Finset.univ.filter fun i => τ < hS.eigenvalues i).card = s)
    {y : Fin s → Fin p → ℝ} {ρ : Fin s → ℝ} (hρ : ∀ k, τ + m ≤ ρ k)
    (hres : ∀ k, (S *ᵥ y k - ρ k • y k) ⬝ᵥ (S *ᵥ y k - ρ k • y k) ≤ δ ^ 2)
    (hgram : ∀ k l, |y k ⬝ᵥ y l - if k = l then 1 else 0| ≤ δ')
    (hsmall : s * (δ' + δ / m) ≤ 1 / 2) (x : Fin p → ℝ) :
    |‖specProj S (Set.Ioi τ) (WithLp.toLp 2 x)‖ ^ 2 - ∑ k, (y k ⬝ᵥ x) ^ 2|
      ≤ 5 * s * (δ' + δ / m) * (x ⬝ᵥ x) := by
  have hnorm : ∀ v : EuclideanSpace ℝ (Fin p), ‖v‖ ^ 2 = WithLp.ofLp v ⬝ᵥ WithLp.ofLp v := by
    intro v
    rw [← real_inner_self_eq_norm_sq, inner_eq_dot]
  have hres' : ∀ k, ‖toOp S (WithLp.toLp 2 (y k)) - ρ k • WithLp.toLp 2 (y k)‖ ≤ δ := by
    intro k
    refine abs_le_of_sq_le_sq ?_ hδ |>.trans' (le_abs_self _)
    rw [hnorm]
    exact hres k
  have hgram' : ∀ k l, |⟪WithLp.toLp 2 (y k), WithLp.toLp 2 (y l)⟫_ℝ - if k = l then 1 else 0|
      ≤ δ' := by
    intro k l
    rw [inner_eq_dot]
    exact hgram k l
  have h := specProj_frame_approx hS hm hδ hδ' hI hρ hres' hgram' hsmall (WithLp.toLp 2 x)
  have hx : ‖(WithLp.toLp 2 x : EuclideanSpace ℝ (Fin p))‖ ^ 2 = x ⬝ᵥ x := by
    rw [hnorm, WithLp.ofLp_toLp]
  rw [hx] at h
  simp only [inner_eq_dot] at h
  exact h

end FrameLemma

/-! ### 6. The split `specProjTop = P_{>τ} + P_edge` (plan 1.3 item 5) -/

section Split

variable {p : ℕ} {S : Matrix (Fin p) (Fin p) ℝ}

/-- When every sorted eigenvalue of index at least `r` is at most `τ` (the count lemma), the
top-`r` eigenvalue set is the disjoint union of the eigenvalues above `τ` and the top-`r`
eigenvalues at most `τ`. -/
theorem specProjTop_eq_specProj_union (hS : S.IsHermitian) {r : ℕ} {τ : ℝ}
    (hτ : ∀ k : Fin (Fintype.card (Fin p)), r ≤ (k : ℕ) → hS.eigenvalues₀ k ≤ τ) :
    specProjTop S hS r = specProj S (Set.Ioi τ ∪ (topEigSet S hS r ∩ Set.Iic τ)) := by
  change specProj S (topEigSet S hS r) = _
  refine specProj_congr_of_iff hS fun i => ?_
  have hi : hS.eigenvalues i = hS.eigenvalues₀ ((eigIdx p).symm i) := by
    rw [← eigenvalues_eigIdx hS, Equiv.apply_symm_apply]
  constructor
  · intro hmem
    by_cases hlt : τ < hS.eigenvalues i
    · exact Or.inl hlt
    · exact Or.inr ⟨hmem, not_lt.mp hlt⟩
  · rintro (hlt | hmem)
    · refine ⟨(eigIdx p).symm i, ?_, hi.symm⟩
      by_contra hge
      have := hτ _ (not_lt.mp hge)
      rw [← hi] at this
      exact absurd hlt (not_lt.mpr this)
    · exact hmem.1

/-- The two eigenvalue sets of the split are disjoint. -/
theorem disjoint_Ioi_edge (hS : S.IsHermitian) (r : ℕ) (τ : ℝ) :
    Disjoint (Set.Ioi τ) (topEigSet S hS r ∩ Set.Iic τ) :=
  Set.disjoint_left.mpr fun _ ht ht' => absurd (Set.mem_Ioi.mp ht) (not_lt.mpr ht'.2)

/-- Plan 1.3 item 5: `specProjTop S r = P_{>τ} + P_edge`, with `P_edge` the projector on the
top-`r` eigenvalues at most `τ`. -/
theorem specProjTop_eq_specProj_add_edge (hS : S.IsHermitian) {r : ℕ} {τ : ℝ}
    (hτ : ∀ k : Fin (Fintype.card (Fin p)), r ≤ (k : ℕ) → hS.eigenvalues₀ k ≤ τ) :
    specProjTop S hS r
      = specProj S (Set.Ioi τ) + specProj S (topEigSet S hS r ∩ Set.Iic τ) := by
  rw [specProjTop_eq_specProj_union hS hτ, specProj_union_of_disjoint hS (disjoint_Ioi_edge hS r τ)]

/-- The squared norms of the split add: `‖P_top x‖² = ‖P_{>τ} x‖² + ‖P_edge x‖²`. -/
theorem norm_sq_specProjTop_split (hS : S.IsHermitian) {r : ℕ} {τ : ℝ}
    (hτ : ∀ k : Fin (Fintype.card (Fin p)), r ≤ (k : ℕ) → hS.eigenvalues₀ k ≤ τ)
    (x : EuclideanSpace ℝ (Fin p)) :
    ‖specProjTop S hS r x‖ ^ 2
      = ‖specProj S (Set.Ioi τ) x‖ ^ 2 + ‖specProj S (topEigSet S hS r ∩ Set.Iic τ) x‖ ^ 2 := by
  rw [specProjTop_eq_specProj_union hS hτ, norm_sq_specProj_union hS (disjoint_Ioi_edge hS r τ)]

end Split

end Frame
end StackedSVD

/-! ### 7. Task U4: the frame approximation with the count derived from the split

Appended 2026-08-30 for task U4 (`notes/archive/plan_subspacelaw.md`, Route P, 1.3 items 3 to 5).
`RankR/RMT/Outliers.lean` consumes this. The count hypothesis `hcount` of
`specProjTop_frame_approx` is derived here from the split `S = W + Q Qᵀ`
(`eigenvalues₀_le_of_split`, at most `r` indices above `τ`) and from the frame itself
(`le_card_filter_of_frame`, at least `r`). -/

namespace StackedSVD
namespace Frame

/-- The count: with the split `S = W + Q Qᵀ`, `lamMax W ≤ τ`, and an `r`-frame of
approximate eigenvectors at points above `τ + mg`, the sorted eigenvalue indices above
`τ` are exactly the top `r`. -/
theorem count_of_split_of_frame {p r : ℕ} {S W : Matrix (Fin p) (Fin p) ℝ}
    {Q : Matrix (Fin p) (Fin r) ℝ} (hW : W.IsHermitian) (hS : S.IsHermitian)
    (hSeq : S = W + Q * Qᵀ) (hrp : r ≤ p) {τ mg δ δ' : ℝ}
    (hτ : lamMax W hW ≤ τ) (hmg : 0 < mg) (hδ' : 0 ≤ δ')
    {y : Fin r → EuclideanSpace ℝ (Fin p)} {ρ : Fin r → ℝ} (hρ : ∀ k, τ + mg ≤ ρ k)
    (hres : ∀ k, ‖toOp S (y k) - ρ k • y k‖ ≤ δ)
    (hgram : ∀ k l, |⟪y k, y l⟫_ℝ - if k = l then 1 else 0| ≤ δ')
    (hsmall : (r : ℝ) * (δ' + (δ / mg) ^ 2) < 1) :
    ∀ k, τ < hS.eigenvalues₀ k ↔ (k : ℕ) < r := by
  classical
  set b : Fin r → Fin p → ℝ := fun k i => ⟪hS.eigenvectorBasis i, y k⟫_ℝ with hb
  have hres' : ∀ k, ∑ i, (hS.eigenvalues i - ρ k) ^ 2 * b k i ^ 2 ≤ δ ^ 2 := by
    intro k
    have h := pow_le_pow_left₀ (norm_nonneg _) (hres k) 2
    rwa [norm_sq_residual hS] at h
  have hgram' : ∀ k l, |b k ⬝ᵥ b l - if k = l then 1 else 0| ≤ δ' := by
    intro k l
    have h : b k ⬝ᵥ b l = ⟪y k, y l⟫_ℝ := by
      rw [inner_eq_sum_coord hS]
      rfl
    rw [h]
    exact hgram k l
  have hlow : r ≤ (Finset.univ.filter fun i => τ < hS.eigenvalues i).card :=
    le_card_filter_of_frame hmg hδ' hρ hres' hgram' hsmall
  have hup : ∀ k : Fin (Fintype.card (Fin p)), r ≤ (k : ℕ) → hS.eigenvalues₀ k ≤ τ :=
    fun k hk => eigenvalues₀_le_of_split hW hS hSeq hτ k hk
  have hsorted : r ≤ (Finset.univ.filter fun k => τ < hS.eigenvalues₀ k).card := by
    rw [← card_filter_eigenvalues hS (fun t => τ < t)]
    exact hlow
  have hsub : (Finset.univ.filter fun k => τ < hS.eigenvalues₀ k)
      ⊆ Finset.univ.filter fun k : Fin (Fintype.card (Fin p)) => (k : ℕ) < r := by
    intro k hk
    rw [Finset.mem_filter] at hk ⊢
    refine ⟨hk.1, ?_⟩
    by_contra hge
    exact absurd hk.2 (not_lt.mpr (hup k (not_lt.mp hge)))
  have hcards : (Finset.univ.filter
      fun k : Fin (Fintype.card (Fin p)) => (k : ℕ) < r).card = r :=
    card_filter_lt (by simpa using hrp)
  have heq := Finset.eq_of_subset_of_card_le hsub (by rw [hcards]; exact hsorted)
  intro k
  constructor
  · intro h
    have hmem : k ∈ Finset.univ.filter fun k => τ < hS.eigenvalues₀ k :=
      Finset.mem_filter.mpr ⟨Finset.mem_univ _, h⟩
    rw [heq] at hmem
    exact (Finset.mem_filter.mp hmem).2
  · intro h
    have hmem : k ∈ Finset.univ.filter
        fun k : Fin (Fintype.card (Fin p)) => (k : ℕ) < r :=
      Finset.mem_filter.mpr ⟨Finset.mem_univ _, h⟩
    rw [← heq] at hmem
    exact (Finset.mem_filter.mp hmem).2

/-- Task U4 assembly of the frame lemma: `specProjTop` against the frame, with the count
derived from the split and the frame instead of supplied. `hδm : δ ≤ mg` turns the frame
smallness into the count smallness. -/
theorem specProjTop_frame_approx_of_split {p r : ℕ} {S W : Matrix (Fin p) (Fin p) ℝ}
    {Q : Matrix (Fin p) (Fin r) ℝ} (hW : W.IsHermitian) (hS : S.IsHermitian)
    (hSeq : S = W + Q * Qᵀ) (hrp : r ≤ p) {τ mg δ δ' : ℝ}
    (hτ : lamMax W hW ≤ τ) (hmg : 0 < mg) (hδ : 0 ≤ δ) (hδ' : 0 ≤ δ')
    {y : Fin r → EuclideanSpace ℝ (Fin p)} {ρ : Fin r → ℝ} (hρ : ∀ k, τ + mg ≤ ρ k)
    (hres : ∀ k, ‖toOp S (y k) - ρ k • y k‖ ≤ δ)
    (hgram : ∀ k l, |⟪y k, y l⟫_ℝ - if k = l then 1 else 0| ≤ δ')
    (hsmall : (r : ℝ) * (δ' + δ / mg) ≤ 1 / 2) (hδm : δ ≤ mg)
    (x : EuclideanSpace ℝ (Fin p)) :
    |‖specProjTop S hS r x‖ ^ 2 - ∑ k, ⟪y k, x⟫_ℝ ^ 2|
      ≤ 5 * r * (δ' + δ / mg) * ‖x‖ ^ 2 := by
  have hdm : δ / mg ≤ 1 := (div_le_one hmg).mpr hδm
  have hdm0 : 0 ≤ δ / mg := div_nonneg hδ hmg.le
  have hr0 : (0 : ℝ) ≤ r := Nat.cast_nonneg r
  have hsmall2 : (r : ℝ) * (δ' + (δ / mg) ^ 2) < 1 := by
    nlinarith [mul_nonneg hr0 (mul_nonneg hdm0 (sub_nonneg.mpr hdm))]
  have hcount := count_of_split_of_frame hW hS hSeq hrp hτ hmg hδ' hρ hres hgram hsmall2
  exact specProjTop_frame_approx hS hrp hcount hmg hδ hδ' hρ hres hgram hsmall x

end Frame
end StackedSVD
