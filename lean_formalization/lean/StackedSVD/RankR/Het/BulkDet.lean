/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.Het.Split
import StackedSVD.RankR.Het.Duality
import StackedSVD.RankR.Het.Simplicity
import StackedSVD.RankR.Het.Deloc
import StackedSVD.RankR.Het.Edge
import StackedSVD.RankR.RMT.EdgeDetR
import StackedSVD.RankR.RMT.EdgeGlueDetR
import StackedSVD.RankR.RMT.DecompR
import StackedSVD.RankR.RMT.EdgeR
import StackedSVD.RankR.RMT.R6R
import StackedSVD.LinAlg.Eigen

/-!
# Stage E8bd: deterministic lemmas for the heteroscedastic bulk step

Track E, stage E8b, deterministic part. Six items of linear algebra that the bulk branch of
`thm:rank_r_stacksvd` (`main_paper.tex:2337`, heteroscedastic, rank `r`) reads at one sample
point. The rank-one mirror is `thm:stacksvd_weighted` (`:463`). No
probability except the two almost-sure statements of section 3, which only reindex the
simplicity facts of stage E5 (`RankR/Het/Simplicity.lean`) and of stage E6
(`RankR/Het/Deloc.lean`).

1. **PSD interlacing** (`eigenvalues₀_le_of_add_posSemidef`). `λ_k(A) ≤ λ_k(A + P)` for every
   positive semidefinite `P`. It generalizes `EdgeDetR.eigenvalues₀_le_of_split`, which is the
   case `P = Q Qᵀ`. The rank-1 mirror is `Frame.eigenvalues₀_le_of_split`.
2. **Residual of a column combination** (`rvecOfR_mulVec`). `rvecOfR Q (Q a) = 0` when
   `Qᵀ Q` is invertible.
3. **Simplicity transfer** between `Xᵀ X` and `X Xᵀ`. `SimpleSpec (Xᵀ X) (rk + 1)` gives
   `SimpleSpec (X Xᵀ) rk` when `rk + 1 ≤ min p q`; the extra index guards against a zero
   eigenvalue in the top block. Two model corollaries give `SimpleSpec` of `X_w X_wᵀ` and of
   `W₀'` almost surely. Rank-1 mirror: `MultiTableModel.topSimple_ae_stackGramW`.
4. **Lower bound of the top-`r` eigenvalues** (`le_eigVal_transpose_mul_self_of_split`).
   `μ₀ ≤ λ_a(Xᵀ X)` for `a < r` when `X Xᵀ = W + Q Qᵀ`, `W` PSD, and `μ₀ Iᵣ ≼ Qᵀ Q`.
5. **Scaling of the test vector** (`normSq_specProj_eq_normSq_mul`) and the corollary of the
   core `EdgeGlueDetR.normSq_specProj_edge_le` at `v = Q a` (`normSq_specProj_edge_le_mulVec`),
   where the residual term of the core is `0`.
6. **Overlap bound through the duality** (`overlapIdx_le_of_eigVal_ge`).
   `overlapIdx X a y ≤ ‖P_a(X Xᵀ) (X y)‖² / μ₀` when `μ₀ ≤ λ_a(Xᵀ X)` and `0 < μ₀`.

Paper labels: items 4 to 6 are the deterministic half of the bulk branch of
`thm:rank_r_stacksvd` (`main_paper.tex:2337`, Appendix E, `sec:rank_r`), whose rank-one
mirror is `thm:stacksvd_weighted` (`:463`). Item 1 carries no paper label. It is
Courant-Fischer monotonicity of the sorted eigenvalues under a positive semidefinite
addition, a generic linear-algebra step, and the paper proves nothing of the kind at this
point: `lem:general_rank_delocalization` (`main_paper.tex:1948`, Appendix D) is the
svdstack delocalization lemma, its proof (`:1957` to `:1974`) has no interlacing step, and
Weyl's inequality enters later, at `:1984`, inside the proof of
`prop:general_rank_unweighted_svdstack`. Rank-1 mirrors in Lean are `SVDStack/Gram.lean`
(gram duality) and `RMT/R6.lean` (the subcritical align bound).

## Conventions

`X : Matrix (Fin p) (Fin q) ℝ` as in `RankR/Het/Duality.lean`: `Xᵀ X` is `q × q` and
`X Xᵀ` is `p × p`, and `min p q` is the number of nonzero singular values.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace HetBulkDet

/-! ### 1. PSD interlacing -/

section Interlace

variable {d : ℕ}

/-- Adding a positive semidefinite matrix never lowers a Rayleigh quotient. It generalizes
`EdgeDetR.rayleighQuotient_le_of_split` (`RankR/RMT/EdgeDetR.lean:123`) from `P = Q Qᵀ`. -/
theorem rayleighQuotient_le_of_add_posSemidef {A P : Matrix (Fin d) (Fin d) ℝ}
    (hP : P.PosSemidef) (x : EuclideanSpace ℝ (Fin d)) :
    LinearMap.rayleighQuotient (toOp A) x ≤ LinearMap.rayleighQuotient (toOp (A + P)) x := by
  have hA : ⟪toOp A x, x⟫_ℝ = (A *ᵥ WithLp.ofLp x) ⬝ᵥ WithLp.ofLp x :=
    Frame.inner_eq_dot _ _
  have hS : ⟪toOp (A + P) x, x⟫_ℝ = ((A + P) *ᵥ WithLp.ofLp x) ⬝ᵥ WithLp.ofLp x :=
    Frame.inner_eq_dot _ _
  have hnn : 0 ≤ (P *ᵥ WithLp.ofLp x) ⬝ᵥ WithLp.ofLp x := by
    have h := hP.dotProduct_mulVec_nonneg (WithLp.ofLp x)
    have hstar : star (WithLp.ofLp x) = WithLp.ofLp x := rfl
    rw [hstar, dotProduct_comm] at h
    exact h
  have hinner : ⟪toOp A x, x⟫_ℝ ≤ ⟪toOp (A + P) x, x⟫_ℝ := by
    rw [hA, hS, Matrix.add_mulVec, add_dotProduct]
    linarith
  unfold LinearMap.rayleighQuotient
  have hre : ∀ a : ℝ, RCLike.re (a : ℝ) = a := fun a => rfl
  simp only [hre]
  have hmul := mul_le_mul_of_nonneg_right hinner (inv_nonneg.mpr (sq_nonneg ‖x‖))
  simpa [div_eq_mul_inv] using hmul

/-- **PSD interlacing, one direction.** `λ_k(A) ≤ λ_k(A + P)` at every sorted index when `P`
is positive semidefinite. It generalizes `EdgeDetR.eigenvalues₀_le_of_split`
(`RankR/RMT/EdgeDetR.lean:145`), the case `P = Q Qᵀ`; the proof is the same Courant-Fischer
count (StatsMLlib `exists_ne_zero_mem_inf_trailingEigenSubspace_of_finrank_eq_succ`). The
statement is false for a merely Hermitian `P`: `P = -A` sends every eigenvalue to `0`. -/
theorem eigenvalues₀_le_of_add_posSemidef {A P : Matrix (Fin d) (Fin d) ℝ}
    (hA : A.IsHermitian) (hP : P.PosSemidef) (k : Fin (Fintype.card (Fin d))) :
    hA.eigenvalues₀ k ≤ (hA.add hP.1).eigenvalues₀ k := by
  have hS : (A + P).IsHermitian := hA.add hP.1
  have hn : Module.finrank ℝ (EuclideanSpace ℝ (Fin d)) = Fintype.card (Fin d) :=
    finrank_euclideanSpace
  have hTA := Frame.isSymmetric_toOp hA
  have hTS := Frame.isSymmetric_toOp hS
  have hkle : (k : ℕ) + 1 ≤ Fintype.card (Fin d) := k.2
  have hLdim : Module.finrank ℝ (hTA.leadingEigenSubspace hn hkle) = (k : ℕ) + 1 :=
    hTA.finrank_leadingEigenSubspace hn hkle
  obtain ⟨x, hxL, hxT, hx0⟩ :=
    hTS.exists_ne_zero_mem_inf_trailingEigenSubspace_of_finrank_eq_succ hn k _ hLdim
  have h1 : hTA.eigenvalues hn k ≤ LinearMap.rayleighQuotient (toOp A) x :=
    hTA.eigenvalues_le_rayleighQuotient_of_mem_leadingEigenSubspace hn k hxL hx0
  have h3 : LinearMap.rayleighQuotient (toOp (A + P)) x ≤ hTS.eigenvalues hn k :=
    hTS.rayleighQuotient_le_eigenvalues_of_mem_trailingEigenSubspace hn k hxT hx0
  have h2 := rayleighQuotient_le_of_add_posSemidef (A := A) hP x
  rw [EdgeDetR.eigenvalues₀_eq_op hA, EdgeDetR.eigenvalues₀_eq_op hS]
  linarith

/-- PSD interlacing on a named sum `S = A + P`, with any proof of `S.IsHermitian`. -/
theorem eigenvalues₀_le_of_eq_add_posSemidef {A P S : Matrix (Fin d) (Fin d) ℝ}
    (hA : A.IsHermitian) (hP : P.PosSemidef) (hS : S.IsHermitian) (hSeq : S = A + P)
    (k : Fin (Fintype.card (Fin d))) :
    hA.eigenvalues₀ k ≤ hS.eigenvalues₀ k := by
  subst hSeq
  exact eigenvalues₀_le_of_add_posSemidef hA hP k

end Interlace

/-! ### 2. The residual of a column combination -/

section Residual

variable {d r : ℕ}

/-- `rvecOfR Q (Q a) = 0`: a vector in the column span of `Q` has no residual, once `Qᵀ Q`
is invertible. Deterministic; `DecompR.aOfR` recovers the coefficient `a` itself. -/
theorem rvecOfR_mulVec {Q : Matrix (Fin d) (Fin r) ℝ} (hQ : IsUnit (Qᵀ * Q).det)
    (a : Fin r → ℝ) : DecompR.rvecOfR Q (Q *ᵥ a) = 0 := by
  have h1 : Qᵀ *ᵥ (Q *ᵥ a) = (Qᵀ * Q) *ᵥ a := Matrix.mulVec_mulVec _ _ _
  have haOf : DecompR.aOfR Q (Q *ᵥ a) = a := by
    rw [DecompR.aOfR, h1, Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hQ,
      Matrix.one_mulVec]
  rw [DecompR.rvecOfR, haOf, sub_self]

end Residual

/-! ### 3. Simplicity transfer between `Xᵀ X` and `X Xᵀ` -/

section Simplicity

variable {p q : ℕ}

/-- `SimpleSpec` at `rk` is `SimpleIdx` at every index below `rk`. The forward direction is
`simpleIdx_of_simpleSpec` (`RankR/Het/Simplicity.lean`); the backward one is definitional. -/
theorem simpleSpec_iff_forall_simpleIdx {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (rk : ℕ) : SimpleSpec A hA rk ↔ ∀ j, j < rk → SimpleIdx A hA j := by
  constructor
  · intro hs j hj
    exact simpleIdx_of_simpleSpec hs hj
  · intro h k l hk hkl
    exact h k hk k l rfl hkl

/-- **The top `rk` eigenvalues of `Xᵀ X` are positive under `SimpleSpec` at `rk + 1`.** The
eigenvalue at index `rk` exists (`rk + 1 ≤ q`) and is nonnegative; the one at `j < rk` is at
least as large and differs from it, so it is positive. Without the extra index a zero
eigenvalue of multiplicity one could sit among the top `rk`. -/
theorem eigVal_pos_of_simpleSpec_succ (X : Matrix (Fin p) (Fin q) ℝ) {rk : ℕ}
    (hs : SimpleSpec (Xᵀ * X) (isHermitian_transpose_mul_self X) (rk + 1))
    (hrk : rk + 1 ≤ min p q) {j : ℕ} (hj : j < rk) :
    0 < eigVal (Xᵀ * X) (isHermitian_transpose_mul_self X) j := by
  have hA := isHermitian_transpose_mul_self X
  have hjq : j < Fintype.card (Fin q) := by rw [Fintype.card_fin]; omega
  have hrq : rk < Fintype.card (Fin q) := by rw [Fintype.card_fin]; omega
  have hnn : 0 ≤ hA.eigenvalues₀ ⟨rk, hrq⟩ := by
    have hpsd : (Xᵀ * X).PosSemidef := by
      simpa using Matrix.posSemidef_conjTranspose_mul_self X
    rw [← eigenvalues_eigIdx]
    exact hpsd.eigenvalues_nonneg _
  have hle : hA.eigenvalues₀ ⟨rk, hrq⟩ ≤ hA.eigenvalues₀ ⟨j, hjq⟩ :=
    hA.eigenvalues₀_antitone (Fin.mk_le_mk.mpr hj.le)
  have hne : hA.eigenvalues₀ ⟨rk, hrq⟩ ≠ hA.eigenvalues₀ ⟨j, hjq⟩ :=
    hs ⟨j, hjq⟩ ⟨rk, hrq⟩ (by simp only []; omega) (by rw [Ne, Fin.ext_iff]; simp only []; omega)
  rw [eigVal_eq _ _ hjq]
  exact lt_of_le_of_lt hnn (lt_of_le_of_ne hle hne)

/-- **Simplicity transfer.** `SimpleSpec (Xᵀ X) (rk + 1)` gives `SimpleSpec (X Xᵀ) rk` when
`rk + 1 ≤ min p q`. Each index `j < rk` has a positive eigenvalue by
`eigVal_pos_of_simpleSpec_succ`, and `Het.simpleIdx_gram_comm` moves `SimpleIdx` across the
duality at such an index. -/
theorem simpleSpec_mul_transpose_self_of_transpose_mul_self (X : Matrix (Fin p) (Fin q) ℝ)
    {rk : ℕ} (hs : SimpleSpec (Xᵀ * X) (isHermitian_transpose_mul_self X) (rk + 1))
    (hrk : rk + 1 ≤ min p q) :
    SimpleSpec (X * Xᵀ) (isHermitian_mul_transpose_self X) rk := by
  rw [simpleSpec_iff_forall_simpleIdx]
  intro j hj
  have hpos := eigVal_pos_of_simpleSpec_succ X hs hrk hj
  exact (Het.simpleIdx_gram_comm X hpos).mp (simpleIdx_of_simpleSpec hs (by omega))

end Simplicity

/-! ### 4. The lower bound of the top eigenvalues -/

section LowerBound

variable {p q rr : ℕ}

/-- A uniform lower bound of the quadratic form bounds every sorted eigenvalue: test the form
at the unit eigenvector `u_k`. -/
theorem le_eigenvalues₀_of_form {A : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian) {μ₀ : ℝ}
    (h : ∀ y : Fin p → ℝ, μ₀ * (y ⬝ᵥ y) ≤ y ⬝ᵥ (A *ᵥ y)) (k : Fin (Fintype.card (Fin p))) :
    μ₀ ≤ hA.eigenvalues₀ k := by
  set i := eigIdx p k with hidef
  have hev : A *ᵥ WithLp.ofLp (hA.eigenvectorBasis i)
      = hA.eigenvalues i • WithLp.ofLp (hA.eigenvectorBasis i) :=
    hA.mulVec_eigenvectorBasis i
  have hnorm : WithLp.ofLp (hA.eigenvectorBasis i) ⬝ᵥ WithLp.ofLp (hA.eigenvectorBasis i)
      = 1 := by
    rw [← Frame.inner_eq_dot, real_inner_self_eq_norm_sq,
      hA.eigenvectorBasis.orthonormal.1 i, one_pow]
  have hk := h (WithLp.ofLp (hA.eigenvectorBasis i))
  rw [hev, dotProduct_smul, smul_eq_mul, hnorm, mul_one, mul_one] at hk
  rw [← eigenvalues_eigIdx hA k]
  exact hk

/-- `μ₀ Iᵣ ≼ Qᵀ Q` as a form bound, from positive semidefiniteness of `Qᵀ Q - μ₀ • 1`. This
is the bridge from `EdgeR.posSemidef_of_close_to_diag` to the `hQQ` hypothesis of the core. -/
theorem form_le_of_posSemidef_sub {G : Matrix (Fin rr) (Fin rr) ℝ} {μ₀ : ℝ}
    (h : (G - μ₀ • (1 : Matrix (Fin rr) (Fin rr) ℝ)).PosSemidef) (y : Fin rr → ℝ) :
    μ₀ * (y ⬝ᵥ y) ≤ y ⬝ᵥ (G *ᵥ y) := by
  have h1 := h.dotProduct_mulVec_nonneg y
  have hstar : star y = y := rfl
  rw [hstar, Matrix.sub_mulVec, dotProduct_sub, Matrix.smul_mulVec, Matrix.one_mulVec,
    dotProduct_smul, smul_eq_mul, sub_nonneg] at h1
  exact h1

/-- `Qᵀ Q - μ₀ I` is PSD when `Qᵀ Q` is within `η` of a diagonal whose entries are at least
`μ₀ + η rr`. The entrywise route of E8b section 2 to the `hQQ` hypothesis; the work is
`EdgeR.posSemidef_of_close_to_diag` (`RankR/RMT/EdgeR.lean:96`) on `G - μ₀ • 1`. -/
theorem posSemidef_sub_smul_one_of_close_to_diag {G : Matrix (Fin rr) (Fin rr) ℝ}
    (hsymm : Gᵀ = G) {L : Fin rr → ℝ} {η μ₀ : ℝ} (hη : 0 ≤ η)
    (hL : ∀ k, μ₀ + η * rr ≤ L k)
    (hclose : ∀ k l, |G k l - (if k = l then L k else 0)| ≤ η) :
    (G - μ₀ • (1 : Matrix (Fin rr) (Fin rr) ℝ)).PosSemidef := by
  have hsymm' : (G - μ₀ • (1 : Matrix (Fin rr) (Fin rr) ℝ))ᵀ = G - μ₀ • 1 := by
    rw [Matrix.transpose_sub, hsymm, Matrix.transpose_smul, Matrix.transpose_one]
  refine EdgeR.posSemidef_of_close_to_diag hsymm' (L := fun k => L k - μ₀) hη
    (fun k => by have := hL k; linarith) fun k l => ?_
  have h := hclose k l
  rw [Matrix.sub_apply, Matrix.smul_apply, Matrix.one_apply, smul_eq_mul]
  by_cases hkl : k = l
  · subst hkl
    simp only [if_true] at h ⊢
    have hre : G k k - μ₀ * 1 - (L k - μ₀) = G k k - L k := by ring
    rw [hre]
    exact h
  · simp only [hkl, if_false] at h ⊢
    simpa using h

/-- The form bound `μ₀ ‖y‖² ≤ yᵀ G y` from entrywise closeness to a diagonal. -/
theorem form_le_of_close_to_diag {G : Matrix (Fin rr) (Fin rr) ℝ} (hsymm : Gᵀ = G)
    {L : Fin rr → ℝ} {η μ₀ : ℝ} (hη : 0 ≤ η) (hL : ∀ k, μ₀ + η * rr ≤ L k)
    (hclose : ∀ k l, |G k l - (if k = l then L k else 0)| ≤ η) (y : Fin rr → ℝ) :
    μ₀ * (y ⬝ᵥ y) ≤ y ⬝ᵥ (G *ᵥ y) :=
  form_le_of_posSemidef_sub (posSemidef_sub_smul_one_of_close_to_diag hsymm hη hL hclose) y

/-- **Item 4, matrix form.** For `X Xᵀ = W + Q Qᵀ` with `W` PSD and `μ₀ Iᵣ ≼ Qᵀ Q`, every
index `a < rr` has `μ₀ ≤ λ_a(Xᵀ X)`. Chain: `λ_a(Xᵀ X) = λ_a(X Xᵀ)` (duality, needs
`a < min p q`), `≥ λ_a(Q Qᵀ)` (PSD interlacing with `P = W`), `= λ_a(Qᵀ Q)` (duality, needs
`a < min p rr`), `≥ μ₀` (the form at the eigenvector). The size hypotheses `rr ≤ p` and
`rr ≤ q` keep every index inside both Gram matrices. -/
theorem le_eigVal_transpose_mul_self_of_split (X : Matrix (Fin p) (Fin q) ℝ)
    {W : Matrix (Fin p) (Fin p) ℝ} {Q : Matrix (Fin p) (Fin rr) ℝ} (hW : W.PosSemidef)
    (hSeq : X * Xᵀ = W + Q * Qᵀ) {μ₀ : ℝ}
    (hQQ : ∀ y : Fin rr → ℝ, μ₀ * (y ⬝ᵥ y) ≤ y ⬝ᵥ ((Qᵀ * Q) *ᵥ y))
    {a : ℕ} (ha : a < rr) (hrp : rr ≤ p) (hrq : rr ≤ q) :
    μ₀ ≤ eigVal (Xᵀ * X) (isHermitian_transpose_mul_self X) a := by
  have hap : a < Fintype.card (Fin p) := by rw [Fintype.card_fin]; omega
  have har : a < Fintype.card (Fin rr) := by rw [Fintype.card_fin]; omega
  have hQ := isHermitian_mul_transpose_self Q
  have hX := isHermitian_mul_transpose_self X
  have h1 : eigVal (Xᵀ * X) (isHermitian_transpose_mul_self X) a
      = eigVal (X * Xᵀ) hX a := Het.eigVal_gram_comm X (by omega)
  have h2 : eigVal (Qᵀ * Q) (isHermitian_transpose_mul_self Q) a
      = eigVal (Q * Qᵀ) hQ a := Het.eigVal_gram_comm Q (by omega)
  have h3 : hQ.eigenvalues₀ ⟨a, hap⟩ ≤ hX.eigenvalues₀ ⟨a, hap⟩ :=
    eigenvalues₀_le_of_eq_add_posSemidef hQ hW hX (by rw [hSeq, add_comm]) ⟨a, hap⟩
  have h4 : μ₀ ≤ eigVal (Qᵀ * Q) (isHermitian_transpose_mul_self Q) a := by
    rw [eigVal_eq _ _ har]
    exact le_eigenvalues₀_of_form _ hQQ ⟨a, har⟩
  rw [h2, eigVal_eq _ _ hap] at h4
  rw [h1, eigVal_eq _ _ hap]
  linarith

/-- Item 4 with the Loewner hypothesis `(Qᵀ Q - μ₀ I).PosSemidef` in place of the form
bound. -/
theorem le_eigVal_transpose_mul_self_of_split' (X : Matrix (Fin p) (Fin q) ℝ)
    {W : Matrix (Fin p) (Fin p) ℝ} {Q : Matrix (Fin p) (Fin rr) ℝ} (hW : W.PosSemidef)
    (hSeq : X * Xᵀ = W + Q * Qᵀ) {μ₀ : ℝ}
    (hQQ : (Qᵀ * Q - μ₀ • (1 : Matrix (Fin rr) (Fin rr) ℝ)).PosSemidef)
    {a : ℕ} (ha : a < rr) (hrp : rr ≤ p) (hrq : rr ≤ q) :
    μ₀ ≤ eigVal (Xᵀ * X) (isHermitian_transpose_mul_self X) a :=
  le_eigVal_transpose_mul_self_of_split X hW hSeq (form_le_of_posSemidef_sub hQQ) ha hrp hrq

/-- Item 4 on `eigenvalues₀` at a `Fin` index below `rr`. -/
theorem le_eigenvalues₀_transpose_mul_self_of_split (X : Matrix (Fin p) (Fin q) ℝ)
    {W : Matrix (Fin p) (Fin p) ℝ} {Q : Matrix (Fin p) (Fin rr) ℝ} (hW : W.PosSemidef)
    (hSeq : X * Xᵀ = W + Q * Qᵀ) {μ₀ : ℝ}
    (hQQ : ∀ y : Fin rr → ℝ, μ₀ * (y ⬝ᵥ y) ≤ y ⬝ᵥ ((Qᵀ * Q) *ᵥ y))
    (k : Fin (Fintype.card (Fin q))) (hk : (k : ℕ) < rr) (hrp : rr ≤ p) (hrq : rr ≤ q) :
    μ₀ ≤ (isHermitian_transpose_mul_self X).eigenvalues₀ k := by
  have h := le_eigVal_transpose_mul_self_of_split X hW hSeq hQQ hk hrp hrq
  rwa [eigVal_eq _ _ k.2] at h

end LowerBound

/-! ### 5. Scaling of the test vector and the core at `v = Q a` -/

section Scaling

variable {p : ℕ}

/-- `‖P v‖² = ‖v‖² ‖P (v / ‖v‖)‖²` for the spectral projector; at `v = 0` both sides are
`0`. -/
theorem normSq_specProj_eq_normSq_mul (S : Matrix (Fin p) (Fin p) ℝ) (E : Set ℝ)
    (v : EuclideanSpace ℝ (Fin p)) :
    ‖specProj S E v‖ ^ 2 = ‖v‖ ^ 2 * ‖specProj S E (‖v‖⁻¹ • v)‖ ^ 2 := by
  by_cases hv : v = 0
  · subst hv
    simp
  · have hn : ‖v‖ ≠ 0 := norm_ne_zero_iff.mpr hv
    rw [map_smul, norm_smul, norm_inv, norm_norm, mul_pow, inv_pow,
      mul_inv_cancel_left₀ (pow_ne_zero 2 hn)]

/-- The same for `specProjIdx`. -/
theorem normSq_specProjIdx_eq_normSq_mul (S : Matrix (Fin p) (Fin p) ℝ) (hS : S.IsHermitian)
    (k : ℕ) (v : EuclideanSpace ℝ (Fin p)) :
    ‖specProjIdx S hS k v‖ ^ 2 = ‖v‖ ^ 2 * ‖specProjIdx S hS k (‖v‖⁻¹ • v)‖ ^ 2 :=
  normSq_specProj_eq_normSq_mul S _ v

/-- **The core of gap G1 at a column combination `v = Q a`.** The residual of `v` along the
columns of `Q` is `0` (`rvecOfR_mulVec`), so `hrest` of
`EdgeGlueDetR.normSq_specProj_edge_le` (`RankR/RMT/EdgeGlueDetR.lean:450`) holds at `ξ = 0`
for the unit vector `v / ‖v‖`, and `normSq_specProj_eq_normSq_mul` restores the scale. The
bound is `‖Q a‖² · 2 rr / (bracket · μQ)`, with the same eight inputs as the core. -/
theorem normSq_specProj_edge_le_mulVec {D rr : ℕ} {W₀ S : Matrix (Fin D) (Fin D) ℝ}
    {Q : Matrix (Fin D) (Fin rr) ℝ}
    (hW₀ : W₀.IsHermitian) (hS : S.IsHermitian) (hSeq : S = W₀ + Q * Qᵀ)
    (hsimpleW : SimpleSpec W₀ hW₀ rr) (hsimpleS : SimpleSpec S hS rr) (hrD : rr ≤ D)
    {z₀ ε₀ τ κ μmin μQ : ℝ}
    (hε₀ : 0 < ε₀) (hedge : lamMax W₀ hW₀ + ε₀ ≤ z₀) (hτ : τ ≤ z₀)
    (hκ : ∑ a ∈ Finset.range rr, ∑ l, ‖specProjIdx W₀ hW₀ a
        (WithLp.toLp 2 fun i => Q i l)‖ ^ 2 ≤ κ)
    (hmin : ∀ y : Fin rr → ℝ, μmin * (y ⬝ᵥ y) ≤ R4.qform2 W₀ z₀ (Q *ᵥ y))
    (hbr : 0 < μmin - (rr : ℝ) * κ / ε₀ ^ 2)
    (hμQ : 0 < μQ) (hQQ : ∀ y : Fin rr → ℝ, μQ * (y ⬝ᵥ y) ≤ y ⬝ᵥ ((Qᵀ * Q) *ᵥ y))
    (a : Fin rr → ℝ) :
    ‖specProj S (topEigSet S hS rr ∩ Set.Iic τ) (WithLp.toLp 2 (Q *ᵥ a))‖ ^ 2
      ≤ ‖(WithLp.toLp 2 (Q *ᵥ a) : EuclideanSpace ℝ (Fin D))‖ ^ 2
          * (2 * ((rr : ℝ) / ((μmin - (rr : ℝ) * κ / ε₀ ^ 2) * μQ))) := by
  set v : EuclideanSpace ℝ (Fin D) := WithLp.toLp 2 (Q *ᵥ a) with hvdef
  by_cases hv0 : v = 0
  · rw [hv0, map_zero, norm_zero]
    simp
  · have hn : ‖v‖ ≠ 0 := norm_ne_zero_iff.mpr hv0
    set u : EuclideanSpace ℝ (Fin D) := ‖v‖⁻¹ • v with hudef
    have hu : ‖u‖ = 1 := by
      rw [hudef, norm_smul, norm_inv, norm_norm, inv_mul_cancel₀ hn]
    have hunit : IsUnit (Qᵀ * Q).det := EdgeGlueDetR.isUnit_det_of_lower_bound Q hμQ hQQ
    have hofLp : WithLp.ofLp u = Q *ᵥ (‖v‖⁻¹ • a) := by
      rw [hudef, WithLp.ofLp_smul, hvdef, WithLp.ofLp_toLp, Matrix.mulVec_smul]
    have hrv : DecompR.rvecOfR Q (WithLp.ofLp u) = 0 := by
      rw [hofLp]
      exact rvecOfR_mulVec hunit _
    have hrest : ∑ k ∈ Finset.range rr,
        ‖specProjIdx S hS k (WithLp.toLp 2 (DecompR.rvecOfR Q (WithLp.ofLp u)))‖ ^ 2 ≤ 0 := by
      rw [hrv]
      simp
    have hcore := EdgeGlueDetR.normSq_specProj_edge_le hW₀ hS hSeq hsimpleW hsimpleS hrD hε₀
      hedge hτ hκ hmin hbr hμQ hQQ hu hrest
    rw [mul_zero, add_zero] at hcore
    rw [normSq_specProj_eq_normSq_mul S _ v]
    exact mul_le_mul_of_nonneg_left hcore (sq_nonneg _)

end Scaling

/-! ### 6. The overlap bound through the duality -/

section Overlap

variable {p q : ℕ}

/-- **Item 6.** When `0 < μ₀ ≤ λ_a(Xᵀ X)` and `a < min p q`, the right overlap at index `a`
is at most the left projection of `X y` divided by `μ₀`. From
`Het.overlapIdx_eq_normSq_specProjIdx_div` (an identity with `λ_a(X Xᵀ)` in the denominator)
and `λ_a(X Xᵀ) = λ_a(Xᵀ X) ≥ μ₀`. The positivity `0 < μ₀` is a hypothesis, not derived. -/
theorem overlapIdx_le_of_eigVal_ge (X : Matrix (Fin p) (Fin q) ℝ) {a : ℕ} (ha : a < min p q)
    {μ₀ : ℝ} (hμ : 0 < μ₀)
    (hle : μ₀ ≤ eigVal (Xᵀ * X) (isHermitian_transpose_mul_self X) a)
    (y : EuclideanSpace ℝ (Fin q)) :
    overlapIdx X a y
      ≤ ‖specProjIdx (X * Xᵀ) (isHermitian_mul_transpose_self X) a
            (WithLp.toLp 2 (X *ᵥ WithLp.ofLp y))‖ ^ 2 / μ₀ := by
  have hpos : 0 < eigVal (Xᵀ * X) (isHermitian_transpose_mul_self X) a := lt_of_lt_of_le hμ hle
  rw [Het.overlapIdx_eq_normSq_specProjIdx_div X hpos y, ← Het.eigVal_gram_comm X ha]
  exact div_le_div_of_nonneg_left (sq_nonneg _) hμ hle

end Overlap

end HetBulkDet

/-! ### 7. The model corollaries -/

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- **`SimpleSpec` of `X_w X_wᵀ` at `r`, almost surely** (item 3, first model corollary).
From `simpleSpec_ae_stackGramW` at `r + 1` on the `d`-side Gram matrix, moved to the row
side by `HetBulkDet.simpleSpec_mul_transpose_self_of_transpose_mul_self`. Needs
`r + 1 ≤ min (∑ n_i) (d N)` and nonzero weights. Rank-1 mirror:
`MultiTableModel.topSimple_ae_stackGramW` (`RMT/Het/Simplicity.lean:203`). -/
theorem simpleSpec_ae_stackXW_mul_transpose [NeZero M] (m : UnalignedModelR μ M n d r rk)
    (w : Fin M → ℝ) (hw : ∀ i, w i ≠ 0) (hG : m.JointGaussianNoise) (N : ℕ)
    (hN : r + 1 ≤ min (∑ i, n i N) (d N)) :
    ∀ᵐ ω ∂(μ N), SimpleSpec (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
      (isHermitian_mul_transpose_self _) r := by
  filter_upwards [m.simpleSpec_ae_stackGramW w hG hw (r + 1) N hN] with ω hω
  exact HetBulkDet.simpleSpec_mul_transpose_self_of_transpose_mul_self (m.stackXW w N ω) hω hN

/-- **Item 3d, `W₀'` half.** `W₀'` has a simple top-`r` spectrum almost surely when
`r ≤ min (∑ n_i) p` and `d N = p + r`. The Gaussian statement is stage E6's
`HetDeloc.simpleSpec_ae_Wsig` (`RankR/Het/Deloc.lean:261`); `W0hetR_eq_Wsig` (`:422`) and the
block law `exists_block_hasLaw_hetR` move it to `μ N` through `ae_of_ae_map`, so no
measurability of the predicate is needed. The case `p = 0` forces `r = 0`, where the
predicate is empty. Rank-1 mirror: stage E5's `MultiTableModel.topSimple_ae_stackGramW` at
the block. -/
theorem simpleSpec_ae_W0hetR [NeZero M] (m : UnalignedModelR μ M n d r rk)
    (w : Fin M → ℝ) (hw : ∀ i, w i ≠ 0) (hG : m.JointGaussianNoise) (N : ℕ) {p : ℕ}
    (hpd : d N = p + r) (hN : r ≤ min (∑ i, n i N) p) :
    ∀ᵐ ω ∂(μ N), SimpleSpec (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) r := by
  rcases Nat.eq_zero_or_pos p with hp0 | hp
  · subst hp0
    have hr : r = 0 := le_antisymm (hN.trans (min_le_right _ _)) (Nat.zero_le _)
    subst hr
    exact Filter.Eventually.of_forall fun ω k l hk _ => absurd hk (Nat.not_lt_zero _)
  obtain ⟨B, hB, hlaw⟩ := m.exists_block_hasLaw_hetR hG N hpd
  have hlawB : HasLaw B (gaussianMatrix (∑ i, n i N) p) (μ N) :=
    measurePreserving_snd.fun_comp_hasLaw hlaw
  have hnN : 0 < ∑ i, n i N := m.stack_row_pos N
  have hdN : 0 < d N := (m.tbl ⟨0, NeZero.pos M⟩).hd N
  have hτ : ∀ j, HetR1.tauOf w (HetR1.blkStack n N) j ≠ 0 := fun j => hw _
  have hg := HetDeloc.simpleSpec_ae_Wsig hnN hp hdN hτ r hN
  rw [← hlawB.map_eq] at hg
  filter_upwards [ae_of_ae_map hlawB.aemeasurable hg] with ω hω
  exact EdgeGlueR.simpleSpec_congr (m.W0hetR_eq_Wsig w N ω (B ω) (hB ω)).symm _ _ hω

/-- **Item 4 in the model.** With `X_w X_wᵀ = W₀' + Q Qᵀ` (`gram_eq_hetR`), `W₀'` PSD
(`posSemidef_W0hetR`), and the form bound `μ₀ Iᵣ ≼ Qᵀ Q`, every index `a < r` has
`μ₀ ≤ λ_a(X_wᵀ X_w)`. The two size hypotheses keep the index inside both Gram matrices. -/
theorem le_eigVal_stackGramW_of_form (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ)
    (N : ℕ) (ω : Ω N) {μ₀ : ℝ}
    (hQQ : ∀ y : Fin r → ℝ, μ₀ * (y ⬝ᵥ y)
      ≤ y ⬝ᵥ (((m.QmatHetR w N ω)ᵀ * m.QmatHetR w N ω) *ᵥ y))
    {a : ℕ} (ha : a < r) (hrn : r ≤ ∑ i, n i N) (hrd : r ≤ d N) :
    μ₀ ≤ eigVal ((m.stackXW w N ω)ᵀ * m.stackXW w N ω)
      (isHermitian_transpose_mul_self _) a :=
  HetBulkDet.le_eigVal_transpose_mul_self_of_split (m.stackXW w N ω) (m.posSemidef_W0hetR w N ω)
    (m.gram_eq_hetR w N ω) hQQ ha hrn hrd

/-- Item 4 for the aligned model, where `r ≤ n_0 N ≤ ∑ n_i` and `r ≤ d N` come from table
`0` (`SpikedModelR.rk_le_n`, `SpikedModelR.rk_le_d`). Stated on `stackGramW`. -/
theorem le_eigVal_stackGramW_of_form_aligned [NeZero M]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N)
    {μ₀ : ℝ}
    (hQQ : ∀ y : Fin r → ℝ, μ₀ * (y ⬝ᵥ y)
      ≤ y ⬝ᵥ (((m.QmatHetR w N ω)ᵀ * m.QmatHetR w N ω) *ᵥ y))
    {a : ℕ} (ha : a < r) :
    μ₀ ≤ eigVal (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) a := by
  have hrn0 : r ≤ n ⟨0, NeZero.pos M⟩ N := (m.tbl ⟨0, NeZero.pos M⟩).rk_le_n N
  have hrn : r ≤ ∑ i, n i N :=
    hrn0.trans (Finset.single_le_sum (f := fun l => n l N) (fun l _ => Nat.zero_le _)
      (Finset.mem_univ (⟨0, NeZero.pos M⟩ : Fin M)))
  have hrd : r ≤ d N := (m.tbl ⟨0, NeZero.pos M⟩).rk_le_d N
  exact m.le_eigVal_stackGramW_of_form w N ω hQQ ha hrn hrd

end UnalignedModelR

end StackedSVD
