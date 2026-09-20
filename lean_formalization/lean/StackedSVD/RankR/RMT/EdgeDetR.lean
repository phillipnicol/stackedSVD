/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.LinAlg.Frame
import StackedSVD.LinAlg.SpecIdx
import StackedSVD.RMT.R4C

/-!
# Task U7c, Part 3: the deterministic bound on `‖Qᵀ u‖²` at an edge eigenvector

Section 4 of `notes/archive/rankr_plan_A.md`, "Part 3, the `col Q` part", and item 3 of section 1.4
of `notes/archive/plan_subspacelaw.md`. The file is pure linear algebra: one real symmetric `W₀`,
one `Q` of `r` columns, `S = W₀ + Q Qᵀ`, and one eigenvector of `S`. No probability and no
asymptotics, so every bound holds for every realization.

The statement. Let `u` be a unit eigenvector of `S` at `lam = λ_k(S)` with `k < r` and
`lam ≤ z₀`, let `y = Qᵀ u`, let `g = z₀ - lamMax W₀ > 0`, let `μmin` be a lower bound of the
form `y ↦ yᵀ Qᵀ G₀(z₀)² Q y` and let `κ` bound `‖Qᵀ u_a(W₀)‖²` over the top `r` eigenvectors
of `W₀`. Then

`(μmin - r κ / g²) ‖y‖² ≤ 1`.

The rank-one mirror is the Cauchy-Schwarz step of `R6.lean` (item R6'), which works with one
column `q` in place of the matrix `Q`.

## Route

1. **Interlacing, one direction** (`eigenvalues₀_le_of_split`). Adding `Q Qᵀ` never lowers a
   sorted eigenvalue: `λ_k(W₀) ≤ λ_k(S)` at every index. The leading eigenspace of `W₀`
   through `k` has dimension `k + 1`, the trailing eigenspace of `S` from `k` has dimension
   `d - k`, so the two meet in a nonzero `x`; the Rayleigh quotient of `x` sits above
   `λ_k(W₀)` on the first, below `λ_k(S)` on the second, and `Q Qᵀ` only raises it. The three
   steps are StatsMLlib's (`LinearAlgebra/Matrix/CourantFischer.lean`); this file supplies the
   bridge `Matrix.IsHermitian.eigenvalues₀ = LinearMap.IsSymmetric.eigenvalues` (it is `rfl`
   on this pin) and the Rayleigh comparison.
2. **The count** (`le_index_of_le_eigenvalues₀`, `card_filter_le_of_simpleSpec`). Under
   `SimpleSpec W₀ hW₀ r` and `λ_k(W₀) ≤ lam` with `k < r`, an index `a` with
   `lam ≤ λ_a(W₀)` has `a ≤ k`. So at most `k + 1 ≤ r` indices survive.
   **Caution:** the plan asked for `Function.Injective hW₀.eigenvalues₀` here. That is false
   when `W₀` is rank deficient (finding U5a, 2026-09-01: the ambient `d × d` Gram of a shorter
   block repeats the eigenvalue `0`), so the hypothesis is `SimpleSpec W₀ hW₀ r`, which
   `SimplicityR.simpleSpec_ae_gaussianMatrix` supplies in both regimes.
3. **The eigen-expansion** (`normSq_transpose_mulVec_le_of_edge`). `(W₀ - lam) u = -Q y`
   gives `(lam - μ_a) ⟨u_a, u⟩ = ⟨u_a, Q y⟩` at every eigenvector `u_a` of `W₀`. Parseval
   turns `‖u‖ = 1` into `∑_a ⟨u_a, u⟩² = 1`; the indices with `μ_a < lam` keep
   `⟨u_a, Q y⟩² / (z₀ - μ_a)²` because `0 < lam - μ_a ≤ z₀ - μ_a`; the full sum over all `a`
   is `qform2 W₀ z₀ (Q y) ≥ μmin ‖y‖²`; and the at most `r` dropped indices each cost at most
   `κ ‖y‖² / g²` by Cauchy-Schwarz.

## Conventions

Eigenvalue indices come in two flavors and both appear. `Fin (Fintype.card (Fin d))` is the
sorted index of `Matrix.IsHermitian.eigenvalues₀`; `Fin d` is the matrix index of
`Matrix.IsHermitian.eigenvalues` and of the sum in `R4.qform2_eq_sum`. `eigIdx d` moves from
the first to the second and `eigenvalues_eigIdx` transports the eigenvalue
(`LinAlg/Eigen.lean`). The hypotheses `hk` and `hκ` are stated in the sorted index, because
that is the index a spike count lives in; the proof works in the matrix index.

Squared norms of coordinate vectors are `‖(WithLp.toLp 2 w : EuclideanSpace ℝ (Fin r))‖ ^ 2`,
not `‖w‖ ^ 2`: the bare `Fin r → ℝ` carries the supremum norm, which would weaken `hκ` by a
factor `r`. `norm_toLp_sq` is the one-line bridge to `w ⬝ᵥ w`.
-/

open Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace EdgeDetR

/-! ### 0. Two bridges

`norm_toLp_sq` reads the Euclidean squared norm of a coordinate vector as a dot product.
`transpose_eigU_mulVec_apply` identifies the coordinate `a` of `R4.eigU`'s change of basis
with the inner product against `Matrix.IsHermitian.eigenvectorBasis a`, which is what lets
`R4.qform2_eq_sum` meet `Frame.norm_sq_eq_sum_coord`. -/

section Bridges

variable {p : ℕ}

/-- The Euclidean squared norm of a coordinate vector is its dot product with itself. -/
theorem norm_toLp_sq (w : Fin p → ℝ) :
    ‖(WithLp.toLp 2 w : EuclideanSpace ℝ (Fin p))‖ ^ 2 = w ⬝ᵥ w := by
  rw [← real_inner_self_eq_norm_sq, Frame.inner_eq_dot]

/-- Over `ℝ`, `IsHermitian` is symmetry of the transpose. -/
theorem transpose_eq_of_isHermitian {A : Matrix (Fin p) (Fin p) ℝ} (h : A.IsHermitian) :
    Aᵀ = A := by
  rwa [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial] at h

/-- The coordinate of `R4.eigU`'s change of basis is the dot product against the eigenvector
of the same matrix index. -/
theorem transpose_eigU_mulVec_apply {W : Matrix (Fin p) (Fin p) ℝ} (hW : W.IsHermitian)
    (w : Fin p → ℝ) (a : Fin p) :
    ((R4.eigU hW)ᵀ *ᵥ w) a = WithLp.ofLp (hW.eigenvectorBasis a) ⬝ᵥ w := by
  simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply, R4.eigU]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Matrix.IsHermitian.eigenvectorUnitary_apply]

end Bridges

/-! ### 1. Interlacing, one direction

`λ_k(W) ≤ λ_k(W + Q Qᵀ)` at every sorted index. The dimension count is StatsMLlib's
`exists_ne_zero_mem_inf_trailingEigenSubspace_of_finrank_eq_succ`; the mirror in this tree is
`Frame.eigenvalues₀_le_of_split`, which bounds the other side under `lamMax W ≤ τ`. -/

section Interlace

variable {d r : ℕ}

/-- `Matrix.IsHermitian.eigenvalues₀` is `LinearMap.IsSymmetric.eigenvalues` of the symmetric
operator. It is `rfl` on this pin, because Mathlib defines the first as the second. -/
theorem eigenvalues₀_eq_op {W : Matrix (Fin d) (Fin d) ℝ} (hW : W.IsHermitian)
    (k : Fin (Fintype.card (Fin d))) :
    hW.eigenvalues₀ k
      = (Frame.isSymmetric_toOp hW).eigenvalues finrank_euclideanSpace k := rfl

/-- Adding `Q Qᵀ` never lowers a Rayleigh quotient. -/
theorem rayleighQuotient_le_of_split {W S : Matrix (Fin d) (Fin d) ℝ}
    {Q : Matrix (Fin d) (Fin r) ℝ} (hSeq : S = W + Q * Qᵀ) (x : EuclideanSpace ℝ (Fin d)) :
    LinearMap.rayleighQuotient (toOp W) x ≤ LinearMap.rayleighQuotient (toOp S) x := by
  have hW : ⟪toOp W x, x⟫_ℝ = (W *ᵥ WithLp.ofLp x) ⬝ᵥ WithLp.ofLp x :=
    Frame.inner_eq_dot _ _
  have hS : ⟪toOp S x, x⟫_ℝ = (S *ᵥ WithLp.ofLp x) ⬝ᵥ WithLp.ofLp x :=
    Frame.inner_eq_dot _ _
  have hextra : (S *ᵥ WithLp.ofLp x) ⬝ᵥ WithLp.ofLp x
      = (W *ᵥ WithLp.ofLp x) ⬝ᵥ WithLp.ofLp x + ∑ j, (Qᵀ *ᵥ WithLp.ofLp x) j ^ 2 := by
    rw [sum_sq_transpose_mulVec Q (WithLp.ofLp x), hSeq, Matrix.add_mulVec, add_dotProduct]
  have hnn : (0 : ℝ) ≤ ∑ j, (Qᵀ *ᵥ WithLp.ofLp x) j ^ 2 :=
    Finset.sum_nonneg fun j _ => sq_nonneg _
  have hinner : ⟪toOp W x, x⟫_ℝ ≤ ⟪toOp S x, x⟫_ℝ := by
    rw [hW, hS, hextra]; linarith
  unfold LinearMap.rayleighQuotient
  have hre : ∀ a : ℝ, RCLike.re (a : ℝ) = a := fun a => rfl
  simp only [hre]
  have hmul := mul_le_mul_of_nonneg_right hinner (inv_nonneg.mpr (sq_nonneg ‖x‖))
  simpa [div_eq_mul_inv] using hmul

/-- **Interlacing, one direction.** For `S = W + Q Qᵀ`, every sorted eigenvalue of `S` is at
least the sorted eigenvalue of `W` at the same index. -/
theorem eigenvalues₀_le_of_split {W S : Matrix (Fin d) (Fin d) ℝ}
    {Q : Matrix (Fin d) (Fin r) ℝ} (hW : W.IsHermitian) (hS : S.IsHermitian)
    (hSeq : S = W + Q * Qᵀ) (k : Fin (Fintype.card (Fin d))) :
    hW.eigenvalues₀ k ≤ hS.eigenvalues₀ k := by
  have hn : Module.finrank ℝ (EuclideanSpace ℝ (Fin d)) = Fintype.card (Fin d) :=
    finrank_euclideanSpace
  have hTW := Frame.isSymmetric_toOp hW
  have hTS := Frame.isSymmetric_toOp hS
  have hkle : (k : ℕ) + 1 ≤ Fintype.card (Fin d) := k.2
  have hLdim : Module.finrank ℝ (hTW.leadingEigenSubspace hn hkle) = (k : ℕ) + 1 :=
    hTW.finrank_leadingEigenSubspace hn hkle
  obtain ⟨x, hxL, hxT, hx0⟩ :=
    hTS.exists_ne_zero_mem_inf_trailingEigenSubspace_of_finrank_eq_succ hn k _ hLdim
  have h1 : hTW.eigenvalues hn k ≤ LinearMap.rayleighQuotient (toOp W) x :=
    hTW.eigenvalues_le_rayleighQuotient_of_mem_leadingEigenSubspace hn k hxL hx0
  have h3 : LinearMap.rayleighQuotient (toOp S) x ≤ hTS.eigenvalues hn k :=
    hTS.rayleighQuotient_le_eigenvalues_of_mem_trailingEigenSubspace hn k hxT hx0
  have h2 := rayleighQuotient_le_of_split hSeq x
  rw [eigenvalues₀_eq_op hW, eigenvalues₀_eq_op hS]
  linarith

end Interlace

/-! ### 2. The count

Under `SimpleSpec W hW r` an eigenvalue at or above `λ_k(W)` with `k < r` has index at most
`k`, so at most `k + 1 ≤ r` indices survive the threshold. -/

section Count

variable {d r : ℕ}

/-- **The index bound.** Under `SimpleSpec W hW r`, an index whose eigenvalue is at or above a
threshold that `λ_k(W)` already meets sits at or before `k`. -/
theorem le_index_of_le_eigenvalues₀ {W : Matrix (Fin d) (Fin d) ℝ} (hW : W.IsHermitian)
    (hsimple : SimpleSpec W hW r) {k a : Fin (Fintype.card (Fin d))} (hk : (k : ℕ) < r)
    {lam : ℝ} (hlam : hW.eigenvalues₀ k ≤ lam) (ha : lam ≤ hW.eigenvalues₀ a) :
    (a : ℕ) ≤ (k : ℕ) := by
  by_contra hcon
  have h : (k : ℕ) < (a : ℕ) := Nat.not_le.mp hcon
  have hka : k ≤ a := Fin.le_def.mpr (le_of_lt h)
  have hanti : hW.eigenvalues₀ a ≤ hW.eigenvalues₀ k := hW.eigenvalues₀_antitone hka
  have hne : k ≠ a := fun hc => by rw [hc] at h; omega
  exact hsimple k a hk hne (le_antisymm hanti (hlam.trans ha))

/-- **The count, sorted index.** At most `r` sorted indices carry an eigenvalue at or above a
threshold that `λ_k(W)` already meets, for some `k < r`. -/
theorem card_filter_le_of_simpleSpec {W : Matrix (Fin d) (Fin d) ℝ} (hW : W.IsHermitian)
    (hsimple : SimpleSpec W hW r) {k : Fin (Fintype.card (Fin d))} (hk : (k : ℕ) < r)
    {lam : ℝ} (hlam : hW.eigenvalues₀ k ≤ lam) :
    ((Finset.univ : Finset (Fin (Fintype.card (Fin d)))).filter
        fun a => lam ≤ hW.eigenvalues₀ a).card ≤ r := by
  classical
  have hle : ((Finset.univ : Finset (Fin (Fintype.card (Fin d)))).filter
      fun a => lam ≤ hW.eigenvalues₀ a).card ≤ (Finset.range r).card := by
    refine Finset.card_le_card_of_injOn (fun a => (a : ℕ)) ?_ ?_
    · intro a ha
      rw [Finset.mem_coe, Finset.mem_filter] at ha
      rw [Finset.mem_coe, Finset.mem_range]
      exact lt_of_le_of_lt (le_index_of_le_eigenvalues₀ hW hsimple hk hlam ha.2) hk
    · intro a _ b _ hab
      exact Fin.ext hab
  simpa using hle

/-- **The count, matrix index.** The same bound in the index the eigenbasis sum uses. -/
theorem card_filter_eigenvalues_le_of_simpleSpec {W : Matrix (Fin d) (Fin d) ℝ}
    (hW : W.IsHermitian) (hsimple : SimpleSpec W hW r)
    {k : Fin (Fintype.card (Fin d))} (hk : (k : ℕ) < r) {lam : ℝ}
    (hlam : hW.eigenvalues₀ k ≤ lam) :
    ((Finset.univ : Finset (Fin d)).filter fun a => lam ≤ hW.eigenvalues a).card ≤ r := by
  classical
  rw [Frame.card_filter_eigenvalues hW (fun s => lam ≤ s)]
  exact card_filter_le_of_simpleSpec hW hsimple hk hlam

end Count

/-! ### 3. The edge bound -/

section Main

variable {d r : ℕ}

/-- **Part 3 of task U7c** (`notes/archive/rankr_plan_A.md` section 4, item 3 of section 1.4 of
`notes/archive/plan_subspacelaw.md`). Deterministic. `u` is a unit eigenvector of
`S = W₀ + Q Qᵀ` at the sorted index `k < r`, its eigenvalue `lam` is at most `z₀`, and `z₀`
sits above the spectrum of `W₀`. Then the mass of `u` on `col Q` obeys

`(μmin - r κ / (z₀ - lamMax W₀)²) ‖Qᵀ u‖² ≤ 1`,

where `μmin` is any lower bound of `y ↦ yᵀ Qᵀ G₀(z₀)² Q y` on the unit sphere and `κ` bounds
`‖Qᵀ u_a(W₀)‖²` over the top `r` eigenvectors of `W₀`.

The rank-one mirror is the Cauchy-Schwarz step of item R6'. The bound is vacuous at a fixed
`z₀`, exactly as there: the caller sends `z₀` down to the bulk edge after the limit in `N`. -/
theorem normSq_transpose_mulVec_le_of_edge
    {W₀ S : Matrix (Fin d) (Fin d) ℝ} (hW₀ : W₀.IsHermitian) (hS : S.IsHermitian)
    {Q : Matrix (Fin d) (Fin r) ℝ} (hSeq : S = W₀ + Q * Qᵀ)
    {z₀ : ℝ} (hz₀ : lamMax W₀ hW₀ < z₀) (hsimple : SimpleSpec W₀ hW₀ r)
    {u : EuclideanSpace ℝ (Fin d)} (hu : ‖u‖ = 1) {lam : ℝ} (heig : toOp S u = lam • u)
    (hk : ∃ k : Fin (Fintype.card (Fin d)), (k : ℕ) < r ∧ hS.eigenvalues₀ k = lam)
    (hlam : lam ≤ z₀) {κ : ℝ}
    (hκ : ∀ a : Fin (Fintype.card (Fin d)), (a : ℕ) < r →
      ‖(WithLp.toLp 2 (Qᵀ *ᵥ WithLp.ofLp (hW₀.eigenvectorBasis (eigIdx d a)))
          : EuclideanSpace ℝ (Fin r))‖ ^ 2 ≤ κ)
    {μmin : ℝ} (hmin : ∀ y : Fin r → ℝ, μmin * (y ⬝ᵥ y) ≤ R4.qform2 W₀ z₀ (Q *ᵥ y)) :
    (μmin - r * κ / (z₀ - lamMax W₀ hW₀) ^ 2)
        * ‖(WithLp.toLp 2 (Qᵀ *ᵥ WithLp.ofLp u) : EuclideanSpace ℝ (Fin r))‖ ^ 2 ≤ 1 := by
  classical
  obtain ⟨k, hkr, hklam⟩ := hk
  set v : Fin d → ℝ := WithLp.ofLp u with hvdef
  set y : Fin r → ℝ := Qᵀ *ᵥ v with hydef
  set g : ℝ := z₀ - lamMax W₀ hW₀ with hgdef
  have hgpos : 0 < g := by rw [hgdef]; linarith
  set mu : Fin d → ℝ := hW₀.eigenvalues with hmudef
  set e : Fin d → (Fin d → ℝ) := fun a => WithLp.ofLp (hW₀.eigenvectorBasis a) with hedef
  set c : Fin d → ℝ := fun a => e a ⬝ᵥ v with hcdef
  set t : Fin d → ℝ := fun a => e a ⬝ᵥ (Q *ᵥ y) with htdef
  set f : Fin d → ℝ := fun a => ((mu a - z₀)⁻¹) ^ 2 * t a ^ 2 with hfdef
  set Y : ℝ := y ⬝ᵥ y with hYdef
  have hYnn : 0 ≤ Y := Frame.dot_self_nonneg y
  -- the target, in dot-product form
  rw [norm_toLp_sq]
  -- κ is nonnegative, because the index `k` itself is admissible
  have hκnn : 0 ≤ κ := le_trans (sq_nonneg _) (hκ k hkr)
  -- (1) the eigenvector equation, in coordinates
  have hSv : S *ᵥ v = lam • v := by
    have h := congrArg WithLp.ofLp heig
    rw [WithLp.ofLp_smul] at h
    exact h
  have hQQ : (Q * Qᵀ) *ᵥ v = Q *ᵥ y := by rw [hydef, ← Matrix.mulVec_mulVec]
  have hWv : W₀ *ᵥ v = lam • v - Q *ᵥ y := by
    refine eq_sub_of_add_eq ?_
    rw [← hQQ, ← Matrix.add_mulVec, ← hSeq]
    exact hSv
  have hWt : W₀ᵀ = W₀ := transpose_eq_of_isHermitian hW₀
  have hcoord : ∀ a, (lam - mu a) * c a = t a := by
    intro a
    have h1 : e a ⬝ᵥ (W₀ *ᵥ v) = mu a * c a := by
      rw [R4.dotProduct_mulVec_comm hWt]
      have hev : W₀ *ᵥ e a = mu a • e a := hW₀.mulVec_eigenvectorBasis a
      rw [hev, dotProduct_smul, smul_eq_mul, hcdef]
      simp only []
      rw [dotProduct_comm]
    have h2 : e a ⬝ᵥ (W₀ *ᵥ v) = lam * c a - t a := by
      rw [hWv, dotProduct_sub, dotProduct_smul, smul_eq_mul, hcdef, htdef]
    rw [h1] at h2
    linarith
  -- (2) Parseval
  have hone : ∑ a, c a ^ 2 = 1 := by
    have h := Frame.norm_sq_eq_sum_coord hW₀ u
    rw [hu, one_pow] at h
    rw [h]
    refine Finset.sum_congr rfl fun a _ => ?_
    simp only [hcdef, hedef, hvdef, Frame.inner_eq_dot]
  -- (3) the squared resolvent form, in the same coordinates
  have hq2 : R4.qform2 W₀ z₀ (Q *ᵥ y) = ∑ a, f a := by
    rw [R4.qform2_eq_sum hW₀ hz₀]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [hfdef, htdef, transpose_eigU_mulVec_apply hW₀]
  -- (4) the interlacing input to the count
  have hklamW : hW₀.eigenvalues₀ k ≤ lam := by
    rw [← hklam]
    exact eigenvalues₀_le_of_split hW₀ hS hSeq k
  -- (5) the retained indices: at most the total mass 1
  have hAle : ∑ a ∈ Finset.univ.filter (fun a : Fin d => ¬ lam ≤ mu a), f a ≤ 1 := by
    have hterm : ∀ a ∈ Finset.univ.filter (fun a : Fin d => ¬ lam ≤ mu a), f a ≤ c a ^ 2 := by
      intro a ha
      rw [Finset.mem_filter] at ha
      have hmlt : mu a < lam := not_le.mp ha.2
      have hppos : 0 < lam - mu a := by linarith
      have hqpos : 0 < z₀ - mu a := by linarith
      have hinv : ((mu a - z₀)⁻¹) ^ 2 = ((z₀ - mu a) ^ 2)⁻¹ := by
        rw [inv_pow]
        congr 1
        ring
      rw [hfdef]
      simp only []
      rw [hinv, inv_mul_le_iff₀ (by positivity : (0 : ℝ) < (z₀ - mu a) ^ 2), ← hcoord a,
        mul_pow]
      have hsq : (lam - mu a) ^ 2 ≤ (z₀ - mu a) ^ 2 := by nlinarith
      exact mul_le_mul_of_nonneg_right hsq (sq_nonneg _)
    calc ∑ a ∈ Finset.univ.filter (fun a : Fin d => ¬ lam ≤ mu a), f a
        ≤ ∑ a ∈ Finset.univ.filter (fun a : Fin d => ¬ lam ≤ mu a), c a ^ 2 :=
          Finset.sum_le_sum hterm
      _ ≤ ∑ a, c a ^ 2 :=
          Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _)
            fun a _ _ => sq_nonneg _
      _ = 1 := hone
  -- (6) the dropped indices: at most `r` of them, each of size `κ Y / g²`
  have hBle : ∑ a ∈ Finset.univ.filter (fun a : Fin d => lam ≤ mu a), f a
      ≤ (r : ℝ) * ((g ^ 2)⁻¹ * (κ * Y)) := by
    have hMnn : (0 : ℝ) ≤ (g ^ 2)⁻¹ * (κ * Y) := by positivity
    have hterm : ∀ a ∈ Finset.univ.filter (fun a : Fin d => lam ≤ mu a),
        f a ≤ (g ^ 2)⁻¹ * (κ * Y) := by
      intro a ha
      rw [Finset.mem_filter] at ha
      have hmge : lam ≤ mu a := ha.2
      have hmle : mu a ≤ lamMax W₀ hW₀ := R4.eigenvalues_le_lamMax hW₀ a
      have hgle : g ≤ z₀ - mu a := by rw [hgdef]; linarith
      have hqpos : 0 < z₀ - mu a := lt_of_lt_of_le hgpos hgle
      -- the resolvent factor
      have hinv : ((mu a - z₀)⁻¹) ^ 2 = ((z₀ - mu a) ^ 2)⁻¹ := by
        rw [inv_pow]; congr 1; ring
      have hinvle : ((mu a - z₀)⁻¹) ^ 2 ≤ (g ^ 2)⁻¹ := by
        rw [hinv]
        exact inv_anti₀ (by positivity) (by nlinarith)
      -- the numerator, by Cauchy-Schwarz and `hκ`
      have hta : t a = (Qᵀ *ᵥ e a) ⬝ᵥ y := by
        rw [htdef]
        simp only []
        rw [Matrix.dotProduct_mulVec, Matrix.mulVec_transpose]
      have hcs : t a ^ 2 ≤ ((Qᵀ *ᵥ e a) ⬝ᵥ (Qᵀ *ᵥ e a)) * Y := by
        rw [hta, hYdef]
        simpa [dotProduct, sq] using
          Finset.sum_mul_sq_le_sq_mul_sq Finset.univ (fun j => (Qᵀ *ᵥ e a) j) fun j => y j
      have hidx : (((eigIdx d).symm a : Fin (Fintype.card (Fin d))) : ℕ) < r := by
        refine lt_of_le_of_lt
          (le_index_of_le_eigenvalues₀ hW₀ hsimple hkr hklamW ?_) hkr
        rw [← eigenvalues_eigIdx hW₀, Equiv.apply_symm_apply]
        exact hmge
      have hkap : (Qᵀ *ᵥ e a) ⬝ᵥ (Qᵀ *ᵥ e a) ≤ κ := by
        have h := hκ ((eigIdx d).symm a) hidx
        rw [norm_toLp_sq, Equiv.apply_symm_apply] at h
        exact h
      have htn : t a ^ 2 ≤ κ * Y := le_trans hcs (mul_le_mul_of_nonneg_right hkap hYnn)
      calc f a = ((mu a - z₀)⁻¹) ^ 2 * t a ^ 2 := rfl
        _ ≤ (g ^ 2)⁻¹ * t a ^ 2 :=
            mul_le_mul_of_nonneg_right hinvle (sq_nonneg _)
        _ ≤ (g ^ 2)⁻¹ * (κ * Y) :=
            mul_le_mul_of_nonneg_left htn (by positivity)
    have hcard : (Finset.univ.filter (fun a : Fin d => lam ≤ mu a)).card ≤ r :=
      card_filter_eigenvalues_le_of_simpleSpec hW₀ hsimple hkr hklamW
    calc ∑ a ∈ Finset.univ.filter (fun a : Fin d => lam ≤ mu a), f a
        ≤ (Finset.univ.filter (fun a : Fin d => lam ≤ mu a)).card • ((g ^ 2)⁻¹ * (κ * Y)) :=
          Finset.sum_le_card_nsmul _ _ _ hterm
      _ = ((Finset.univ.filter (fun a : Fin d => lam ≤ mu a)).card : ℝ)
            * ((g ^ 2)⁻¹ * (κ * Y)) := nsmul_eq_mul _ _
      _ ≤ (r : ℝ) * ((g ^ 2)⁻¹ * (κ * Y)) :=
          mul_le_mul_of_nonneg_right (by exact_mod_cast hcard) hMnn
  -- (7) assemble
  have hsum : ∑ a, f a
      = ∑ a ∈ Finset.univ.filter (fun a : Fin d => lam ≤ mu a), f a
        + ∑ a ∈ Finset.univ.filter (fun a : Fin d => ¬ lam ≤ mu a), f a :=
    (Finset.sum_filter_add_sum_filter_not _ _ _).symm
  have hlow : μmin * Y ≤ ∑ a, f a := by
    rw [← hq2, hYdef]
    exact hmin y
  have hfinal : μmin * Y ≤ (r : ℝ) * ((g ^ 2)⁻¹ * (κ * Y)) + 1 := by
    rw [hsum] at hlow
    linarith
  have hrew : (r : ℝ) * ((g ^ 2)⁻¹ * (κ * Y)) = ((r : ℝ) * κ / g ^ 2) * Y := by
    field_simp
  rw [hrew] at hfinal
  rw [← hYdef]
  have hexp : (μmin - (r : ℝ) * κ / g ^ 2) * Y
      = μmin * Y - ((r : ℝ) * κ / g ^ 2) * Y := by ring
  rw [hexp]
  linarith

end Main

end EdgeDetR

end StackedSVD
