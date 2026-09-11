/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.LinAlg.SpecProjPerturb

/-!
# Ky Fan's trace inequalities

Task T5a of `notes/RANK_R_PLAN.md`. Neither Mathlib nor StatsMLlib has any trace inequality
for symmetric matrices (`grep -ri "ky.?fan"` is empty in both pinned trees, 2026-08-30), so
this file builds one.

Let `M` be a real symmetric `p × p` matrix with sorted eigenvalues
`λ_0 ≥ λ_1 ≥ ... ≥ λ_{p-1}` (`Matrix.IsHermitian.eigenvalues₀`, which is antitone). For every
`Y : Matrix (Fin p) (Fin r) ℝ` with orthonormal columns (`Yᵀ * Y = 1`):

* `kyFan_max`: `tr (Yᵀ M Y) ≤ ∑_{k < r} λ_k`, attained by the top `r` eigenvectors
  (`kyFan_max_attained`, `kyFan_max_isGreatest`).
* `kyFan_min`: `∑_{k < r} λ_{p-1-k} ≤ tr (Yᵀ M Y)`, attained by the bottom `r` eigenvectors
  (`kyFan_min_attained`, `kyFan_min_isLeast`).

The bottom sum is written `∑ k : Fin r, hM.eigenvalues₀ (Fin.castLE hr k).rev`. `Fin.rev`
sends the index `k` to `p - 1 - k`, so the `k = 0` term is the smallest eigenvalue. That is the
paper's `∑_{ℓ=1}^{r} λ_{p+1-ℓ}` in `thm:gen_rank_weight_svdstak`, at `ℓ = k + 1`.

## Route

Write `M = U diag(λ) Uᵀ` in the eigenbasis and put `Z = Uᵀ Y`. Then

  `tr (Yᵀ M Y) = ∑_a λ_a w_a`,  `w_a = ‖row a of Z‖²`  (`trace_eq_sum_frameWeight`),

and `Z Zᵀ` is an orthogonal projector of rank `r`, so `0 ≤ w_a ≤ 1` and `∑_a w_a = r`
(`frameWeight_nonneg`, `frameWeight_le_one`, `sum_frameWeight`). A weighted sum with such
weights is largest when the weight `1` sits on the `r` largest eigenvalues. That last step is
the standalone rearrangement lemma `sum_mul_le_sum_of_threshold` and its twin
`sum_le_sum_mul_of_threshold`; both are stated with an explicit threshold `t` that splits the
selected index set from its complement, which avoids any exchange or induction argument.

## Index convention

Every statement takes `hr : r ≤ Fintype.card (Fin p)`, the index type of `eigenvalues₀`. A
caller who holds `hrp : r ≤ p` supplies `by simpa using hrp`, since `Fintype.card_fin` is a
simp lemma.

## Numeric check

`p = 8`, `r = 3`, 20000 draws of a GOE matrix and a random orthonormal frame (seed 20260830):
0 violations of either bound, weights inside `[4.5e-5, 0.9952]`, `|∑_a w_a - r| ≤ 5.3e-15`,
attainment error `≤ 2.9e-14`. On `M = diag(3, 1, -2)` with `r = 2` the empirical supremum is
3.99988 (exact 4) and the infimum is -0.99920 (exact -1). The congruence form of section 7 was
checked the same way at `p = 6`, `r = 2` with `A` positive definite and `D` positive
semidefinite: 0 violations, attainment error `≤ 6.3e-15`.
-/

open scoped Matrix

namespace StackedSVD

/-! ### 1. The rearrangement lemma on weights

`lam` is a family of reals, `S` a finite index set, and `t` a threshold with `lam ≥ t` on `S`
and `lam ≤ t` off `S`. If `w` takes values in `[0, 1]` and `∑ w = |S|`, then `∑ lam * w` is at
most `∑_S lam`. Nothing here knows about matrices or about sorting. -/

section Rearrangement

variable {ι : Type*} [Fintype ι]

/-- Ky Fan's rearrangement step, maximum form. -/
theorem sum_mul_le_sum_of_threshold (lam w : ι → ℝ) (S : Finset ι) (t : ℝ)
    (hin : ∀ a ∈ S, t ≤ lam a) (hout : ∀ a ∉ S, lam a ≤ t)
    (hw0 : ∀ a, 0 ≤ w a) (hw1 : ∀ a, w a ≤ 1) (hsum : ∑ a, w a = S.card) :
    ∑ a, lam a * w a ≤ ∑ a ∈ S, lam a := by
  classical
  have hsplit : ∑ a, lam a * w a = (∑ a ∈ S, lam a * w a) + ∑ a ∈ Sᶜ, lam a * w a :=
    (Finset.sum_add_sum_compl S _).symm
  have hwtot : (∑ a ∈ S, t * w a) + ∑ a ∈ Sᶜ, t * w a = t * S.card := by
    rw [Finset.sum_add_sum_compl S (fun a => t * w a), ← Finset.mul_sum, hsum]
  have hS : ∑ a ∈ S, lam a * w a ≤ (∑ a ∈ S, lam a) - t * S.card + ∑ a ∈ S, t * w a := by
    have hstep : ∑ a ∈ S, lam a * w a ≤ ∑ a ∈ S, (lam a + (t * w a - t)) := by
      refine Finset.sum_le_sum fun a ha => ?_
      nlinarith [mul_nonneg (sub_nonneg.2 (hin a ha)) (sub_nonneg.2 (hw1 a))]
    refine hstep.trans_eq ?_
    rw [Finset.sum_add_distrib, Finset.sum_sub_distrib, Finset.sum_const, nsmul_eq_mul]
    ring
  have hC : ∑ a ∈ Sᶜ, lam a * w a ≤ ∑ a ∈ Sᶜ, t * w a :=
    Finset.sum_le_sum fun a ha =>
      mul_le_mul_of_nonneg_right (hout a (Finset.mem_compl.mp ha)) (hw0 a)
  rw [hsplit]
  linarith

/-- Ky Fan's rearrangement step, minimum form. -/
theorem sum_le_sum_mul_of_threshold (lam w : ι → ℝ) (S : Finset ι) (t : ℝ)
    (hin : ∀ a ∈ S, lam a ≤ t) (hout : ∀ a ∉ S, t ≤ lam a)
    (hw0 : ∀ a, 0 ≤ w a) (hw1 : ∀ a, w a ≤ 1) (hsum : ∑ a, w a = S.card) :
    ∑ a ∈ S, lam a ≤ ∑ a, lam a * w a := by
  classical
  have hsplit : ∑ a, lam a * w a = (∑ a ∈ S, lam a * w a) + ∑ a ∈ Sᶜ, lam a * w a :=
    (Finset.sum_add_sum_compl S _).symm
  have hwtot : (∑ a ∈ S, t * w a) + ∑ a ∈ Sᶜ, t * w a = t * S.card := by
    rw [Finset.sum_add_sum_compl S (fun a => t * w a), ← Finset.mul_sum, hsum]
  have hS : (∑ a ∈ S, lam a) - t * S.card + ∑ a ∈ S, t * w a ≤ ∑ a ∈ S, lam a * w a := by
    have hstep : ∑ a ∈ S, (lam a + (t * w a - t)) ≤ ∑ a ∈ S, lam a * w a := by
      refine Finset.sum_le_sum fun a ha => ?_
      nlinarith [mul_nonneg (sub_nonneg.2 (hin a ha)) (sub_nonneg.2 (hw1 a))]
    refine le_trans (le_of_eq ?_) hstep
    rw [Finset.sum_add_distrib, Finset.sum_sub_distrib, Finset.sum_const, nsmul_eq_mul]
    ring
  have hC : ∑ a ∈ Sᶜ, t * w a ≤ ∑ a ∈ Sᶜ, lam a * w a :=
    Finset.sum_le_sum fun a ha =>
      mul_le_mul_of_nonneg_right (hout a (Finset.mem_compl.mp ha)) (hw0 a)
  rw [hsplit]
  linarith

end Rearrangement

/-! ### 2. The weight vector of an orthonormal frame -/

section FrameWeight

variable {p r : ℕ}

/-- The squared norm of row `a` of `Z`. For `Z` with orthonormal columns this is the `a`-th
diagonal entry of the orthogonal projector `Z Zᵀ`, hence a number in `[0, 1]`, and the weights
sum to `r`. -/
def frameWeight (Z : Matrix (Fin p) (Fin r) ℝ) (a : Fin p) : ℝ := ∑ j, Z a j ^ 2

theorem frameWeight_nonneg (Z : Matrix (Fin p) (Fin r) ℝ) (a : Fin p) : 0 ≤ frameWeight Z a :=
  Finset.sum_nonneg fun _ _ => sq_nonneg _

/-- `frameWeight Z a` is the `(a, a)` entry of `Z Zᵀ`. -/
theorem frameWeight_eq_mul_transpose (Z : Matrix (Fin p) (Fin r) ℝ) (a : Fin p) :
    frameWeight Z a = (Z * Zᵀ) a a := by
  simp [frameWeight, Matrix.mul_apply, sq]

/-- The weights of an orthonormal frame sum to the number of columns. -/
theorem sum_frameWeight {Z : Matrix (Fin p) (Fin r) ℝ} (hZ : Zᵀ * Z = 1) :
    ∑ a, frameWeight Z a = r := by
  have hdiag : ∀ j : Fin r, ∑ a, Z a j ^ 2 = 1 := by
    intro j
    have h := congrArg (fun N : Matrix (Fin r) (Fin r) ℝ => N j j) hZ
    simpa [Matrix.mul_apply, Matrix.one_apply, sq] using h
  calc ∑ a, frameWeight Z a = ∑ j, ∑ a, Z a j ^ 2 := Finset.sum_comm
    _ = r := by simp [hdiag]

/-- Each weight of an orthonormal frame is at most `1`. -/
theorem frameWeight_le_one {Z : Matrix (Fin p) (Fin r) ℝ} (hZ : Zᵀ * Z = 1) (a : Fin p) :
    frameWeight Z a ≤ 1 := by
  set P : Matrix (Fin p) (Fin p) ℝ := Z * Zᵀ with hP
  have hidem : P * P = P := by
    rw [hP, Matrix.mul_assoc, ← Matrix.mul_assoc Zᵀ Z Zᵀ, hZ, Matrix.one_mul]
  have hsymm : ∀ i j, P j i = P i j := by
    intro i j
    simp [hP, Matrix.mul_apply, mul_comm]
  have hentry : P a a = ∑ b, P a b ^ 2 := by
    have h := congrArg (fun N : Matrix (Fin p) (Fin p) ℝ => N a a) hidem
    simp only [Matrix.mul_apply] at h
    rw [← h]
    exact Finset.sum_congr rfl fun b _ => by rw [hsymm a b]; ring
  have hge : P a a ^ 2 ≤ P a a := by
    conv_rhs => rw [hentry]
    exact Finset.single_le_sum (f := fun b => P a b ^ 2) (fun b _ => sq_nonneg _)
      (Finset.mem_univ a)
  rw [frameWeight_eq_mul_transpose, ← hP]
  nlinarith [hge]

end FrameWeight

/-! ### 3. Diagonalization in the eigenbasis

Three lemmas repeated from `RMT/R4.lean`, which this file cannot import (`R4.lean` sits above
`LinAlg/` in the import graph). They are `private`, so no name clashes with `R4.lean`. -/

section Eigenbasis

variable {p : ℕ} {M : Matrix (Fin p) (Fin p) ℝ}

/-- The orthogonal matrix whose columns are the eigenvectors of `M`. -/
private noncomputable def eigW (hM : M.IsHermitian) : Matrix (Fin p) (Fin p) ℝ :=
  (hM.eigenvectorUnitary : Matrix (Fin p) (Fin p) ℝ)

private theorem transpose_eigW_mul (hM : M.IsHermitian) : (eigW hM)ᵀ * eigW hM = 1 := by
  have h := hM.eigenvectorUnitary.2
  rw [Matrix.mem_unitaryGroup_iff'] at h
  rw [eigW, ← Matrix.conjTranspose_eq_transpose_of_trivial (α := ℝ)]
  exact h

private theorem eigW_mul_transpose (hM : M.IsHermitian) : eigW hM * (eigW hM)ᵀ = 1 := by
  have h := hM.eigenvectorUnitary.2
  rw [Matrix.mem_unitaryGroup_iff] at h
  rw [eigW, ← Matrix.conjTranspose_eq_transpose_of_trivial (α := ℝ)]
  exact h

private theorem eigW_conj (hM : M.IsHermitian) :
    eigW hM * Matrix.diagonal hM.eigenvalues * (eigW hM)ᵀ = M := by
  conv_rhs => rw [hM.spectral_theorem]
  simp [eigW, Unitary.conjStarAlgAut_apply, Matrix.star_eq_conjTranspose,
    Matrix.conjTranspose_eq_transpose_of_trivial]

/-- `W` cancels `Wᵀ` on the left of any matrix of matching size. -/
private theorem eigW_cancel_left {q : ℕ} (hM : M.IsHermitian) (X : Matrix (Fin p) (Fin q) ℝ) :
    eigW hM * ((eigW hM)ᵀ * X) = X := by
  rw [← Matrix.mul_assoc, eigW_mul_transpose hM, Matrix.one_mul]

/-- `Wᵀ` cancels `W` on the left of any matrix of matching size. -/
private theorem transpose_eigW_cancel_left {q : ℕ} (hM : M.IsHermitian)
    (X : Matrix (Fin p) (Fin q) ℝ) : (eigW hM)ᵀ * (eigW hM * X) = X := by
  rw [← Matrix.mul_assoc, transpose_eigW_mul hM, Matrix.one_mul]

/-- `Wᵀ M W` is the diagonal matrix of the eigenvalues. -/
private theorem transpose_eigW_mul_mul (hM : M.IsHermitian) :
    (eigW hM)ᵀ * M * eigW hM = Matrix.diagonal hM.eigenvalues := by
  have h : (eigW hM)ᵀ * (eigW hM * Matrix.diagonal hM.eigenvalues * (eigW hM)ᵀ) * eigW hM
      = (eigW hM)ᵀ * M * eigW hM := by rw [eigW_conj hM]
  rw [← h]
  simp only [Matrix.mul_assoc]
  rw [transpose_eigW_mul hM, Matrix.mul_one, transpose_eigW_cancel_left hM]

variable {r : ℕ}

/-- The rotated frame `Z = Wᵀ Y` is orthonormal when `Y` is. -/
private theorem transpose_rot_mul_rot (hM : M.IsHermitian) {Y : Matrix (Fin p) (Fin r) ℝ}
    (hY : Yᵀ * Y = 1) : ((eigW hM)ᵀ * Y)ᵀ * ((eigW hM)ᵀ * Y) = 1 := by
  rw [Matrix.transpose_mul, Matrix.transpose_transpose, Matrix.mul_assoc,
    ← Matrix.mul_assoc (eigW hM) ((eigW hM)ᵀ) Y, eigW_mul_transpose hM, Matrix.one_mul, hY]

/-- The trace in the eigenbasis: `tr (Yᵀ M Y) = ∑_a λ_a w_a` with `w = frameWeight (Wᵀ Y)`. -/
theorem trace_eq_sum_frameWeight (hM : M.IsHermitian) (Y : Matrix (Fin p) (Fin r) ℝ) :
    Matrix.trace (Yᵀ * M * Y) = ∑ a, hM.eigenvalues a * frameWeight ((eigW hM)ᵀ * Y) a := by
  have hrw : Yᵀ * M * Y
      = ((eigW hM)ᵀ * Y)ᵀ * Matrix.diagonal hM.eigenvalues * ((eigW hM)ᵀ * Y) := by
    rw [← transpose_eigW_mul_mul hM]
    simp only [Matrix.transpose_mul, Matrix.transpose_transpose, Matrix.mul_assoc]
    rw [eigW_cancel_left hM, eigW_cancel_left hM]
  set Z := (eigW hM)ᵀ * Y with hZ
  rw [hrw]
  have hentry : ∀ j : Fin r,
      (Zᵀ * Matrix.diagonal hM.eigenvalues * Z) j j = ∑ a, hM.eigenvalues a * Z a j ^ 2 := by
    intro j
    rw [Matrix.mul_apply]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [Matrix.mul_diagonal, Matrix.transpose_apply]
    ring
  simp only [Matrix.trace, Matrix.diag_apply, hentry]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun a _ => by rw [frameWeight, Finset.mul_sum]

end Eigenbasis

/-! ### 4. The two index sets: the top `r` and the bottom `r` sorted positions -/

section IndexSets

/-- `k ↦ (Fin.castLE hr k).rev`, the bottom `r` sorted positions in decreasing order of
eigenvalue index: `k = 0` gives the last position `n - 1`. -/
private def revEmb {n r : ℕ} (hr : r ≤ n) : Fin r ↪ Fin n :=
  ⟨fun k => (Fin.castLE hr k).rev, fun _ _ h => Fin.castLE_inj.mp (Fin.rev_injective h)⟩

private theorem mem_map_castLEEmb {n r : ℕ} (hr : r ≤ n) (k : Fin n) :
    k ∈ Finset.map (Fin.castLEEmb hr) Finset.univ ↔ (k : ℕ) < r := by
  simp only [Finset.mem_map, Finset.mem_univ, true_and, Fin.coe_castLEEmb]
  constructor
  · rintro ⟨j, rfl⟩
    simp
  · intro h
    exact ⟨⟨(k : ℕ), h⟩, by simp⟩

private theorem mem_map_revEmb {n r : ℕ} (hr : r ≤ n) (k : Fin n) :
    k ∈ Finset.map (revEmb hr) Finset.univ ↔ n - r ≤ (k : ℕ) := by
  have hk : (k : ℕ) < n := k.isLt
  simp only [Finset.mem_map, Finset.mem_univ, true_and]
  constructor
  · rintro ⟨j, rfl⟩
    have hj : (j : ℕ) < r := j.isLt
    change n - r ≤ (((Fin.castLE hr j).rev : Fin n) : ℕ)
    rw [Fin.val_rev, Fin.val_castLE]
    omega
  · intro h
    have hjlt : n - 1 - (k : ℕ) < r := by omega
    refine ⟨⟨n - 1 - (k : ℕ), hjlt⟩, ?_⟩
    have hjv : ((⟨n - 1 - (k : ℕ), hjlt⟩ : Fin r) : ℕ) = n - 1 - (k : ℕ) := rfl
    change ((Fin.castLE hr ⟨n - 1 - (k : ℕ), hjlt⟩).rev : Fin n) = k
    apply Fin.ext
    rw [Fin.val_rev, Fin.val_castLE, hjv]
    omega

private theorem card_map_castLEEmb {n r : ℕ} (hr : r ≤ n) :
    (Finset.map (Fin.castLEEmb hr) (Finset.univ : Finset (Fin r))).card = r := by
  simp

private theorem card_map_revEmb {n r : ℕ} (hr : r ≤ n) :
    (Finset.map (revEmb hr) (Finset.univ : Finset (Fin r))).card = r := by
  simp

end IndexSets

/-! ### 5. Ky Fan's trace inequalities -/

section KyFan

variable {p r : ℕ} {M : Matrix (Fin p) (Fin p) ℝ}

/-- The eigenvalue list indexed by `Fin p` is the sorted list read through the canonical
equivalence `Fin (Fintype.card (Fin p)) ≃ Fin p`. -/
private theorem eigenvalues_equiv (hM : M.IsHermitian)
    (k : Fin (Fintype.card (Fin p))) :
    hM.eigenvalues ((Fintype.equivOfCardEq (Fintype.card_fin _)) k) = hM.eigenvalues₀ k := by
  simp only [Matrix.IsHermitian.eigenvalues, Equiv.symm_apply_apply]

/-- **Ky Fan's maximum principle.** The trace of `M` compressed to any orthonormal `r`-frame is
at most the sum of the `r` largest eigenvalues of `M`. -/
theorem kyFan_max (hM : M.IsHermitian) (hr : r ≤ Fintype.card (Fin p))
    {Y : Matrix (Fin p) (Fin r) ℝ} (hY : Yᵀ * Y = 1) :
    Matrix.trace (Yᵀ * M * Y) ≤ ∑ k : Fin r, hM.eigenvalues₀ (Fin.castLE hr k) := by
  classical
  rcases Nat.eq_zero_or_pos r with hr0 | hr0
  · subst hr0
    simp [Matrix.trace]
  set n := Fintype.card (Fin p) with hn
  set e : Fin n ≃ Fin p := Fintype.equivOfCardEq (Fintype.card_fin _) with he
  set Z := (eigW hM)ᵀ * Y with hZ
  have hZo : Zᵀ * Z = 1 := transpose_rot_mul_rot hM hY
  set w : Fin p → ℝ := frameWeight Z with hw
  set S : Finset (Fin n) := Finset.map (Fin.castLEEmb hr) Finset.univ with hS
  set t : ℝ := hM.eigenvalues₀ ⟨r - 1, by omega⟩ with ht
  have hcard : S.card = r := card_map_castLEEmb hr
  have hsum : ∑ k : Fin n, w (e k) = (S.card : ℝ) := by
    rw [hcard, Equiv.sum_comp e w]
    exact_mod_cast sum_frameWeight hZo
  have hin : ∀ a ∈ S, t ≤ hM.eigenvalues₀ a := by
    intro a ha
    have hlt : (a : ℕ) < r := (mem_map_castLEEmb hr a).mp ha
    exact hM.eigenvalues₀_antitone (Fin.le_def.mpr (by simp; omega))
  have hout : ∀ a ∉ S, hM.eigenvalues₀ a ≤ t := by
    intro a ha
    have hge : r ≤ (a : ℕ) := by
      by_contra hcon
      exact ha ((mem_map_castLEEmb hr a).mpr (by omega))
    exact hM.eigenvalues₀_antitone (Fin.le_def.mpr (by simp; omega))
  have hkey := sum_mul_le_sum_of_threshold (fun k : Fin n => hM.eigenvalues₀ k)
    (fun k : Fin n => w (e k)) S t hin hout
    (fun k => frameWeight_nonneg Z (e k)) (fun k => frameWeight_le_one hZo (e k)) hsum
  have hlhs : ∑ k : Fin n, hM.eigenvalues₀ k * w (e k)
      = Matrix.trace (Yᵀ * M * Y) := by
    rw [trace_eq_sum_frameWeight hM Y, ← hZ, ← hw,
      ← Equiv.sum_comp e (fun a => hM.eigenvalues a * w a)]
    exact Finset.sum_congr rfl fun k _ => by rw [eigenvalues_equiv hM k]
  have hrhs : ∑ a ∈ S, hM.eigenvalues₀ a = ∑ k : Fin r, hM.eigenvalues₀ (Fin.castLE hr k) := by
    rw [hS, Finset.sum_map]
    rfl
  rw [← hlhs, ← hrhs]
  exact hkey

/-- **Ky Fan's minimum principle.** The trace of `M` compressed to any orthonormal `r`-frame is
at least the sum of the `r` smallest eigenvalues of `M`, written in the paper's index form
`∑_{ℓ=1}^{r} λ_{p+1-ℓ}`. -/
theorem kyFan_min (hM : M.IsHermitian) (hr : r ≤ Fintype.card (Fin p))
    {Y : Matrix (Fin p) (Fin r) ℝ} (hY : Yᵀ * Y = 1) :
    ∑ k : Fin r, hM.eigenvalues₀ (Fin.castLE hr k).rev ≤ Matrix.trace (Yᵀ * M * Y) := by
  classical
  rcases Nat.eq_zero_or_pos r with hr0 | hr0
  · subst hr0
    simp [Matrix.trace]
  set n := Fintype.card (Fin p) with hn
  set e : Fin n ≃ Fin p := Fintype.equivOfCardEq (Fintype.card_fin _) with he
  set Z := (eigW hM)ᵀ * Y with hZ
  have hZo : Zᵀ * Z = 1 := transpose_rot_mul_rot hM hY
  set w : Fin p → ℝ := frameWeight Z with hw
  set S : Finset (Fin n) := Finset.map (revEmb hr) Finset.univ with hS
  set t : ℝ := hM.eigenvalues₀ ⟨n - r, by omega⟩ with ht
  have hcard : S.card = r := card_map_revEmb hr
  have hsum : ∑ k : Fin n, w (e k) = (S.card : ℝ) := by
    rw [hcard, Equiv.sum_comp e w]
    exact_mod_cast sum_frameWeight hZo
  have hin : ∀ a ∈ S, hM.eigenvalues₀ a ≤ t := by
    intro a ha
    have hge : n - r ≤ (a : ℕ) := (mem_map_revEmb hr a).mp ha
    exact hM.eigenvalues₀_antitone (Fin.le_def.mpr (by simp; omega))
  have hout : ∀ a ∉ S, t ≤ hM.eigenvalues₀ a := by
    intro a ha
    have hlt : (a : ℕ) < n - r := by
      by_contra hcon
      exact ha ((mem_map_revEmb hr a).mpr (by omega))
    exact hM.eigenvalues₀_antitone (Fin.le_def.mpr (by simp; omega))
  have hkey := sum_le_sum_mul_of_threshold (fun k : Fin n => hM.eigenvalues₀ k)
    (fun k : Fin n => w (e k)) S t hin hout
    (fun k => frameWeight_nonneg Z (e k)) (fun k => frameWeight_le_one hZo (e k)) hsum
  have hrhs : ∑ k : Fin n, hM.eigenvalues₀ k * w (e k)
      = Matrix.trace (Yᵀ * M * Y) := by
    rw [trace_eq_sum_frameWeight hM Y, ← hZ, ← hw,
      ← Equiv.sum_comp e (fun a => hM.eigenvalues a * w a)]
    exact Finset.sum_congr rfl fun k _ => by rw [eigenvalues_equiv hM k]
  have hlhs : ∑ a ∈ S, hM.eigenvalues₀ a
      = ∑ k : Fin r, hM.eigenvalues₀ (Fin.castLE hr k).rev := by
    rw [hS, Finset.sum_map]
    rfl
  rw [← hrhs, ← hlhs]
  exact hkey

end KyFan

/-! ### 6. Attainment -/

section Attainment

variable {p r : ℕ} {M : Matrix (Fin p) (Fin p) ℝ}

/-- Any injective family of `r` eigenvectors is an orthonormal frame, and the compressed trace
is the sum of the corresponding eigenvalues. -/
theorem exists_frame_trace_eq (hM : M.IsHermitian) (σ : Fin r → Fin p)
    (hσ : Function.Injective σ) :
    ∃ Y : Matrix (Fin p) (Fin r) ℝ, Yᵀ * Y = 1 ∧
      Matrix.trace (Yᵀ * M * Y) = ∑ j, hM.eigenvalues (σ j) := by
  classical
  set W := eigW hM with hW
  refine ⟨Matrix.of fun i j => W i (σ j), ?_, ?_⟩
  · ext j k
    have hcomp : ((Matrix.of fun i j => W i (σ j))ᵀ * Matrix.of fun i j => W i (σ j)) j k
        = (Wᵀ * W) (σ j) (σ k) := by
      simp [Matrix.mul_apply, Matrix.transpose_apply]
    rw [hcomp, hW, transpose_eigW_mul hM]
    by_cases h : j = k
    · subst h
      simp
    · rw [Matrix.one_apply_ne (fun hh => h (hσ hh)), Matrix.one_apply_ne h]
  · set Y : Matrix (Fin p) (Fin r) ℝ := Matrix.of fun i j => W i (σ j) with hY
    have hZ : ∀ (a : Fin p) (j : Fin r), (Wᵀ * Y) a j = if a = σ j then (1 : ℝ) else 0 := by
      intro a j
      have hcomp : (Wᵀ * Y) a j = (Wᵀ * W) a (σ j) := by
        simp [hY, Matrix.mul_apply, Matrix.transpose_apply]
      rw [hcomp, hW, transpose_eigW_mul hM, Matrix.one_apply]
    rw [trace_eq_sum_frameWeight hM Y, ← hW]
    have hwt : ∀ a : Fin p, frameWeight (Wᵀ * Y) a = ∑ j, (if a = σ j then (1 : ℝ) else 0) := by
      intro a
      refine Finset.sum_congr rfl fun j _ => ?_
      rw [hZ a j]
      by_cases h : a = σ j <;> simp [h]
    calc ∑ a, hM.eigenvalues a * frameWeight (Wᵀ * Y) a
        = ∑ a, ∑ j, hM.eigenvalues a * (if a = σ j then (1 : ℝ) else 0) := by
          exact Finset.sum_congr rfl fun a _ => by rw [hwt a, Finset.mul_sum]
      _ = ∑ j, ∑ a, hM.eigenvalues a * (if a = σ j then (1 : ℝ) else 0) := Finset.sum_comm
      _ = ∑ j, hM.eigenvalues (σ j) := by
          refine Finset.sum_congr rfl fun j _ => ?_
          simp

/-- The maximum of `kyFan_max` is attained, by the top `r` eigenvectors. -/
theorem kyFan_max_attained (hM : M.IsHermitian) (hr : r ≤ Fintype.card (Fin p)) :
    ∃ Y : Matrix (Fin p) (Fin r) ℝ, Yᵀ * Y = 1 ∧
      Matrix.trace (Yᵀ * M * Y) = ∑ k : Fin r, hM.eigenvalues₀ (Fin.castLE hr k) := by
  classical
  set e : Fin (Fintype.card (Fin p)) ≃ Fin p := Fintype.equivOfCardEq (Fintype.card_fin _) with he
  obtain ⟨Y, hY, htr⟩ := exists_frame_trace_eq hM (fun k : Fin r => e (Fin.castLE hr k))
    (fun a b h => Fin.castLE_inj.mp (e.injective h))
  exact ⟨Y, hY, by rw [htr]; exact Finset.sum_congr rfl fun k _ => eigenvalues_equiv hM _⟩

/-- The minimum of `kyFan_min` is attained, by the bottom `r` eigenvectors. -/
theorem kyFan_min_attained (hM : M.IsHermitian) (hr : r ≤ Fintype.card (Fin p)) :
    ∃ Y : Matrix (Fin p) (Fin r) ℝ, Yᵀ * Y = 1 ∧
      Matrix.trace (Yᵀ * M * Y) = ∑ k : Fin r, hM.eigenvalues₀ (Fin.castLE hr k).rev := by
  classical
  set e : Fin (Fintype.card (Fin p)) ≃ Fin p := Fintype.equivOfCardEq (Fintype.card_fin _) with he
  obtain ⟨Y, hY, htr⟩ := exists_frame_trace_eq hM (fun k : Fin r => e (Fin.castLE hr k).rev)
    (fun a b h => Fin.castLE_inj.mp (Fin.rev_injective (e.injective h)))
  exact ⟨Y, hY, by rw [htr]; exact Finset.sum_congr rfl fun k _ => eigenvalues_equiv hM _⟩

/-- The set of compressed traces over all orthonormal `r`-frames. -/
def traceSet (M : Matrix (Fin p) (Fin p) ℝ) (r : ℕ) : Set ℝ :=
  {t | ∃ Y : Matrix (Fin p) (Fin r) ℝ, Yᵀ * Y = 1 ∧ Matrix.trace (Yᵀ * M * Y) = t}

/-- Ky Fan, maximum form: the sum of the top `r` eigenvalues is the greatest compressed trace. -/
theorem kyFan_max_isGreatest (hM : M.IsHermitian) (hr : r ≤ Fintype.card (Fin p)) :
    IsGreatest (traceSet M r) (∑ k : Fin r, hM.eigenvalues₀ (Fin.castLE hr k)) := by
  refine ⟨?_, ?_⟩
  · obtain ⟨Y, hY, htr⟩ := kyFan_max_attained hM hr
    exact ⟨Y, hY, htr⟩
  · rintro s ⟨Y, hY, rfl⟩
    exact kyFan_max hM hr hY

/-- Ky Fan, minimum form: the sum of the bottom `r` eigenvalues is the least compressed
trace. -/
theorem kyFan_min_isLeast (hM : M.IsHermitian) (hr : r ≤ Fintype.card (Fin p)) :
    IsLeast (traceSet M r) (∑ k : Fin r, hM.eigenvalues₀ (Fin.castLE hr k).rev) := by
  refine ⟨?_, ?_⟩
  · obtain ⟨Y, hY, htr⟩ := kyFan_min_attained hM hr
    exact ⟨Y, hY, htr⟩
  · rintro s ⟨Y, hY, rfl⟩
    exact kyFan_min hM hr hY

end Attainment

/-! ### 7. The congruence form, for `thm:gen_rank_weight_svdstak` (task T5b)

The paper needs the extremes of `tr(Xᵀ D X)` under the constraint `Xᵀ A X = I_r`, with `A`
positive definite. The substitution `Y = S X`, where `S` is a symmetric invertible square root
of `A`, turns the constraint into `Yᵀ Y = I_r` and the objective into `tr(Yᵀ (S⁻¹ D S⁻¹) Y)`,
so both extremes come from section 5.

`S` is a hypothesis, not `CFC.sqrt A`. `Matrix.PosDef.sqrt` and `Matrix.PosSemidef.sqrt` no
longer exist on this Mathlib pin (only `Matrix.PosSemidef.inv_sqrt` survives, in
`Analysis/Matrix/Order.lean:131`, stated for `CFC.sqrt`). Taking `S` as an argument keeps this
file free of the `MatrixOrder` scoped instances, and the consumer supplies
`S := CFC.sqrt A` with `CFC.sq_sqrt` for `S * S = A`, or any other symmetric root. -/

section Congruence

variable {p r : ℕ}

/-- A symmetric matrix has a symmetric inverse. -/
private theorem transpose_inv_eq {S : Matrix (Fin p) (Fin p) ℝ} (hS : Sᵀ = S) : S⁻¹ᵀ = S⁻¹ := by
  rw [Matrix.transpose_nonsing_inv, hS]

/-- `S⁻¹ D S⁻¹` is symmetric when `D` and `S` are. -/
theorem isHermitian_inv_mul_mul_inv {D S : Matrix (Fin p) (Fin p) ℝ} (hD : D.IsHermitian)
    (hS : Sᵀ = S) : (S⁻¹ * D * S⁻¹).IsHermitian := by
  have hDt : Dᵀ = D := by
    rwa [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial] at hD
  rw [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial]
  rw [Matrix.transpose_mul, Matrix.transpose_mul, transpose_inv_eq hS, hDt, Matrix.mul_assoc]

/-- The two feasible sets carry the same traces: `X ↦ S X` is a bijection from the frames with
`Xᵀ A X = 1` to the frames with `Yᵀ Y = 1`, and it matches `tr(Xᵀ D X)` with
`tr(Yᵀ (S⁻¹ D S⁻¹) Y)`. -/
theorem traceSet_congr_eq {A D S : Matrix (Fin p) (Fin p) ℝ} (hS : Sᵀ = S) (hSA : S * S = A)
    (hSu : IsUnit S.det) :
    {t | ∃ X : Matrix (Fin p) (Fin r) ℝ, Xᵀ * A * X = 1 ∧ Matrix.trace (Xᵀ * D * X) = t}
      = traceSet (S⁻¹ * D * S⁻¹) r := by
  have hSi : S⁻¹ * S = 1 := Matrix.nonsing_inv_mul S hSu
  have hiS : S * S⁻¹ = 1 := Matrix.mul_nonsing_inv S hSu
  have c1 : ∀ X : Matrix (Fin p) (Fin r) ℝ, S⁻¹ * (S * X) = X := fun X => by
    rw [← Matrix.mul_assoc, hSi, Matrix.one_mul]
  have c2 : ∀ X : Matrix (Fin p) (Fin r) ℝ, S * (S⁻¹ * X) = X := fun X => by
    rw [← Matrix.mul_assoc, hiS, Matrix.one_mul]
  ext t
  constructor
  · rintro ⟨X, hX, rfl⟩
    refine ⟨S * X, ?_, ?_⟩
    · rw [Matrix.transpose_mul, hS, ← hX]
      simp only [Matrix.mul_assoc, ← hSA]
    · congr 1
      rw [Matrix.transpose_mul, hS]
      simp only [Matrix.mul_assoc]
      rw [c2, c1]
  · rintro ⟨Y, hY, rfl⟩
    refine ⟨S⁻¹ * Y, ?_, ?_⟩
    · rw [Matrix.transpose_mul, transpose_inv_eq hS, ← hY, ← hSA]
      simp only [Matrix.mul_assoc]
      rw [c2, c1]
    · congr 1
      rw [Matrix.transpose_mul, transpose_inv_eq hS]
      simp only [Matrix.mul_assoc]

/-- **Ky Fan, congruence form, maximum.** With `A = S * S` positive definite,
`max { tr(Xᵀ D X) : Xᵀ A X = I_r }` is the sum of the top `r` eigenvalues of `S⁻¹ D S⁻¹`. -/
theorem kyFan_max_congr {A D S : Matrix (Fin p) (Fin p) ℝ} (hS : Sᵀ = S) (hSA : S * S = A)
    (hSu : IsUnit S.det) (hM : (S⁻¹ * D * S⁻¹).IsHermitian) (hr : r ≤ Fintype.card (Fin p)) :
    IsGreatest
      {t | ∃ X : Matrix (Fin p) (Fin r) ℝ, Xᵀ * A * X = 1 ∧ Matrix.trace (Xᵀ * D * X) = t}
      (∑ k : Fin r, hM.eigenvalues₀ (Fin.castLE hr k)) := by
  rw [traceSet_congr_eq hS hSA hSu]
  exact kyFan_max_isGreatest hM hr

/-- **Ky Fan, congruence form, minimum.** This is the half that `thm:gen_rank_weight_svdstak`
uses: `min { tr(Xᵀ D X) : Xᵀ A X = I_r } = ∑_{ℓ=1}^{r} λ_{p+1-ℓ}(S⁻¹ D S⁻¹)`. -/
theorem kyFan_min_congr {A D S : Matrix (Fin p) (Fin p) ℝ} (hS : Sᵀ = S) (hSA : S * S = A)
    (hSu : IsUnit S.det) (hM : (S⁻¹ * D * S⁻¹).IsHermitian) (hr : r ≤ Fintype.card (Fin p)) :
    IsLeast
      {t | ∃ X : Matrix (Fin p) (Fin r) ℝ, Xᵀ * A * X = 1 ∧ Matrix.trace (Xᵀ * D * X) = t}
      (∑ k : Fin r, hM.eigenvalues₀ (Fin.castLE hr k).rev) := by
  rw [traceSet_congr_eq hS hSA hSu]
  exact kyFan_min_isLeast hM hr

end Congruence

end StackedSVD
