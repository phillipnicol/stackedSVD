/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.Stack
import StackedSVD.RankR.RMT.ShiftR
import StackedSVD.RankR.General
import StackedSVD.LinAlg.Eigen

/-!
# Task C1: one rank-`rk` table as a `RankRStack`

Task C1 of `notes/archive/rankr_plan_C.md` section 2. The Gaussian chain of Section 7 (`EdgeR`,
`AlignOutG`, `Forms`, `EdgeGlueR`) reads a `RankRStack` (`RankR/RMT/Stack.lean`), and the
target of Track C, `SpikedModelR.TableLawR` (`RankR/General.lean`), reads a `SpikedModelR`.
This file is the bridge between the two. The rank-1 mirror is `UnalignedModel.toStack`
(`RankR/RMT/Split.lean:657`), which embeds the stacked multi-table model in the same
interface.

## Content

1. `SpikedModelR.toStack`: `X = U Θ Vᵀ + E` is the stack with signal factor `A = U Θ`,
   diagonal core `C = AᵀA = Θ²` and the same noise. Every derived object is a `rfl` bridge
   (`toStack_gram`, `toStack_core`, `overlapIdx_toStack`).
2. `SpikedModelR.rk_le_n`, `rk_le_d`: `rk = rank 1 = rank (UᵀU) ≤ rank U ≤ n N`, and the same
   for `V`. Item C1.2; `simpleSpec_ae_affine` (task C7) takes both.
3. The spectrum of the diagonal core. Mathlib has no formula for
   `IsHermitian.eigenvalues (diagonal g)`, so `coreEig` and `spikeVec` do not read the model
   index `k` but Mathlib's own eigen-index. `eigenvalues₀_core` reads the sorted list off the
   characteristic polynomial (`eigenvalues₀_eq_of_charpoly` with `Matrix.charpoly_diagonal`;
   the `Antitone` side condition is `hθanti` with `hθnn`), and `coreIdx`, `corePerm` carry
   the two index changes. `exists_perm_coreEig` (item C1.3) is the export:
   `coreEig (σ k) = θ k ^ 2` and `spikeVec (σ k) N = ± col N k`.
4. `normSq_specProj_spikeVec_eq` (item C1.4): the sign cancels, because `specProj` is linear
   and the statement squares the norm.
5. `SpikedModelR.shift`, the mirror of `RankRStack.shift` (`RankR/RMT/ShiftR.lean:55`) and of
   `SpikedModel.shift` (`RMT/TailShift.lean:48`), with the `rfl` bridges,
   `shift_Regime`, `shift_GaussianNoise`, `toStack_shift` and the transfer
   `tendstoInProb_overlapIdx_of_shift`. Task C5.2 removes the side condition
   `∀ N, n N = rk + pp N` with it, through `RankRStack.exists_shift_lt`.

Nothing here is random and nothing here is a paper theorem: this is the plumbing tasks C5,
C6, C7 and C9 stand on.

The duplicate `rk_le_d` is removed; the one copy lives in `RankR/General.lean`
(F28, 2026-09-08).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace SpikedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n d : ℕ → ℕ} {rk : ℕ}

/-! ### 1. The embedding -/

/-- **Item C1.1.** One rank-`rk` table is a `RankRStack` with signal factor `A = U Θ` and
diagonal core `C = Θ²`. The mirror at `r_i = 1` is `UnalignedModel.toStack`. -/
noncomputable def toStack (m : SpikedModelR μ n d rk) : RankRStack μ n d rk where
  V := m.V
  A := fun N => m.U N * Matrix.diagonal m.θ
  core := Matrix.diagonal fun k => m.θ k ^ 2
  Zu := m.Z
  E := m.E
  X := m.X
  hV := m.hV
  hcore := by
    intro N
    rw [Matrix.transpose_mul, Matrix.diagonal_transpose, Matrix.mul_assoc,
      ← Matrix.mul_assoc (m.U N)ᵀ, m.hU N, Matrix.one_mul, Matrix.diagonal_mul_diagonal]
    simp [sq]
  hE := fun _ _ => rfl
  hX := fun _ _ => rfl
  hd := m.hd
  hZmeas := m.hZ

section Bridges

variable (m : SpikedModelR μ n d rk)

@[simp] theorem toStack_V (N : ℕ) : m.toStack.V N = m.V N := rfl

@[simp] theorem toStack_core : m.toStack.core = Matrix.diagonal fun k => m.θ k ^ 2 := rfl

@[simp] theorem toStack_X (N : ℕ) (ω : Ω N) : m.toStack.X N ω = m.X N ω := rfl

@[simp] theorem toStack_E (N : ℕ) (ω : Ω N) : m.toStack.E N ω = m.E N ω := rfl

@[simp] theorem toStack_gram (N : ℕ) (ω : Ω N) :
    m.toStack.gram N ω = (m.X N ω)ᵀ * m.X N ω := rfl

/-- The noise law is the same predicate on both sides. -/
theorem gaussianNoise_toStack (hG : m.GaussianNoise) : m.toStack.GaussianNoise := hG

/-- The overlap of `TableLawR` is the index projector of the stack Gram matrix. -/
theorem overlapIdx_toStack (N : ℕ) (ω : Ω N) (k : ℕ) (w : EuclideanSpace ℝ (Fin (d N))) :
    overlapIdx (m.X N ω) k w
      = ‖specProjIdx (m.toStack.gram N ω) (m.toStack.isHermitian_gram N ω) k w‖ ^ 2 := rfl

end Bridges

/-! ### 2. The two rank facts -/

/-- **Item C1.2.** `rk ≤ n N`: the columns of `U N` are orthonormal, so
`rk = rank 1 = rank (UᵀU) ≤ rank U ≤ n N`. -/
theorem rk_le_n (m : SpikedModelR μ n d rk) (N : ℕ) : rk ≤ n N := by
  have h1 : ((m.U N)ᵀ * m.U N).rank ≤ (m.U N).rank := Matrix.rank_mul_le_right _ _
  rw [m.hU N, Matrix.rank_one, Fintype.card_fin] at h1
  have h2 : (m.U N).rank ≤ n N := by
    simpa using Matrix.rank_le_card_height (m.U N)
  exact h1.trans h2

/-! ### 3. The spectrum of the diagonal core -/

/-- The index change of `Matrix.IsHermitian.eigenvalues`, at the core of the embedding: the
model index of the eigenvalue Mathlib stores at the index `j`. -/
noncomputable def coreIdx (rk : ℕ) : Equiv.Perm (Fin rk) :=
  (eigIdx rk).symm.trans (finCongr (Fintype.card_fin rk))

/-- The inverse of `coreIdx`: the eigen-index of the `k`-th spike. This is the permutation
`exists_perm_coreEig` returns. -/
noncomputable def corePerm (rk : ℕ) : Equiv.Perm (Fin rk) := (coreIdx rk).symm

/-- The squares of the signal strengths are pairwise distinct: `hθanti` with `hθnn`. -/
theorem thetaSq_injective (m : SpikedModelR μ n d rk) :
    Function.Injective fun k => m.θ k ^ 2 := by
  intro a b hab
  by_contra hne
  rcases lt_or_gt_of_ne hne with h | h
  · have hlt := m.hθanti h
    nlinarith [m.hθnn a, m.hθnn b]
  · have hlt := m.hθanti h
    nlinarith [m.hθnn a, m.hθnn b]

/-- **The sorted spectrum of the core.** `C = Θ²` is diagonal with a strictly decreasing
nonnegative diagonal, so its characteristic polynomial names the sorted list
(`eigenvalues₀_eq_of_charpoly`, `LinAlg/Eigen.lean`). -/
theorem eigenvalues₀_core (m : SpikedModelR μ n d rk) :
    m.toStack.isHermitian_core.eigenvalues₀
      = fun q => m.θ (finCongr (Fintype.card_fin rk) q) ^ 2 := by
  refine eigenvalues₀_eq_of_charpoly _ ?_ ?_
  · intro a b hab
    have hle : (finCongr (Fintype.card_fin rk)) a ≤ (finCongr (Fintype.card_fin rk)) b :=
      Fin.le_def.mpr (by simpa using Fin.le_def.mp hab)
    have hθ := m.hθanti.antitone hle
    nlinarith [m.hθnn ((finCongr (Fintype.card_fin rk)) a),
      m.hθnn ((finCongr (Fintype.card_fin rk)) b)]
  · rw [m.toStack_core, Matrix.charpoly_diagonal]
    exact (Fintype.prod_equiv (finCongr (Fintype.card_fin rk)) _ _ fun q => rfl).symm

/-- `coreEig j` is the square of the signal strength at the model index `coreIdx rk j`. -/
theorem coreEig_eq (m : SpikedModelR μ n d rk) (j : Fin rk) :
    m.toStack.coreEig j = m.θ (coreIdx rk j) ^ 2 := by
  have h := eigenvalues_eigIdx m.toStack.isHermitian_core ((eigIdx rk).symm j)
  rw [Equiv.apply_symm_apply] at h
  rw [RankRStack.coreEig, h, m.eigenvalues₀_core]
  rfl

/-- The eigenvector of the diagonal core at the eigen-index `j` is the standard basis vector
of the model index `coreIdx rk j`, up to sign. -/
theorem exists_sign_eigenvectorBasis_core (m : SpikedModelR μ n d rk) (j : Fin rk) :
    ∃ t : ℝ, (t = 1 ∨ t = -1) ∧
      WithLp.ofLp (m.toStack.isHermitian_core.eigenvectorBasis j)
        = fun l => if l = coreIdx rk j then t else 0 := by
  set H := m.toStack.isHermitian_core with hH
  set x : Fin rk → ℝ := WithLp.ofLp (H.eigenvectorBasis j) with hx
  set i : Fin rk := coreIdx rk j with hi
  -- The eigenvalue equation, entrywise on the diagonal core.
  have hmul : m.toStack.core *ᵥ x = H.eigenvalues j • x := H.mulVec_eigenvectorBasis j
  have hentry : ∀ l : Fin rk, m.θ l ^ 2 * x l = m.θ i ^ 2 * x l := by
    intro l
    have hlhs : (m.toStack.core *ᵥ x) l = m.θ l ^ 2 * x l := by
      rw [m.toStack_core]
      exact Matrix.mulVec_diagonal _ _ _
    have hval : H.eigenvalues j = m.θ i ^ 2 := m.coreEig_eq j
    have h := congrFun hmul l
    rw [hlhs, Pi.smul_apply, smul_eq_mul, hval] at h
    exact h
  -- Off the index `i` the coordinate vanishes, because the squares are distinct.
  have hzero : ∀ l : Fin rk, l ≠ i → x l = 0 := by
    intro l hl
    have h := hentry l
    have hne : m.θ l ^ 2 - m.θ i ^ 2 ≠ 0 := by
      intro h0
      exact hl (m.thetaSq_injective (by linarith [sub_eq_zero.mp h0]))
    have : (m.θ l ^ 2 - m.θ i ^ 2) * x l = 0 := by linarith
    rcases mul_eq_zero.mp this with h1 | h1
    · exact absurd h1 hne
    · exact h1
  -- The vector is a unit vector, so the surviving coordinate is `± 1`.
  have hone : ‖(H.eigenvectorBasis j : EuclideanSpace ℝ (Fin rk))‖ = 1 :=
    H.eigenvectorBasis.orthonormal.1 j
  have hsum : ∑ l, x l ^ 2 = 1 := by
    have h2 := EuclideanSpace.norm_eq (H.eigenvectorBasis j)
    rw [hone] at h2
    have h3 : ∑ l, ‖x l‖ ^ 2 = 1 := Real.sqrt_eq_one.mp h2.symm
    simpa [Real.norm_eq_abs, sq_abs] using h3
  have hxi : x i ^ 2 = 1 := by
    rw [← hsum]
    refine (Finset.sum_eq_single (f := fun l => x l ^ 2) i (fun l _ hl => ?_)
      (fun h => absurd (Finset.mem_univ i) h)).symm
    rw [hzero l hl]
    ring
  refine ⟨x i, ?_, ?_⟩
  · have : (x i - 1) * (x i + 1) = 0 := by nlinarith [hxi]
    rcases mul_eq_zero.mp this with h1 | h1
    · exact Or.inl (by linarith)
    · exact Or.inr (by linarith)
  · funext l
    by_cases hl : l = i
    · rw [hl, if_pos rfl]
    · rw [if_neg hl, hzero l hl]

/-- The spike direction of the stack at the eigen-index `j` is the column of `V` at the model
index `coreIdx rk j`, up to sign. -/
theorem spikeVec_eq_or (m : SpikedModelR μ n d rk) (j : Fin rk) (N : ℕ) :
    m.toStack.spikeVec j N = m.col N (coreIdx rk j)
      ∨ m.toStack.spikeVec j N = -m.col N (coreIdx rk j) := by
  obtain ⟨t, ht, hb⟩ := m.exists_sign_eigenvectorBasis_core j
  have hmv : (m.V N).mulVec (WithLp.ofLp (m.toStack.isHermitian_core.eigenvectorBasis j))
      = t • fun l => m.V N l (coreIdx rk j) := by
    funext l
    rw [Matrix.mulVec, dotProduct]
    simp only [hb, mul_ite, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, if_true,
      Pi.smul_apply, smul_eq_mul]
    ring
  have hspike : m.toStack.spikeVec j N = t • m.col N (coreIdx rk j) := by
    rw [RankRStack.spikeVec, m.toStack_V, hmv]
    rfl
  rcases ht with rfl | rfl
  · exact Or.inl (by rw [hspike, one_smul])
  · exact Or.inr (by rw [hspike, neg_smul, one_smul])

/-- **Item C1.3.** The permutation that matches the model index `k` with the eigen-index
`σ k` of the core: the eigenvalue is `θ k ^ 2` and the spike direction is the `k`-th column of
`V`, up to a sign that `normSq_specProj_spikeVec_eq` cancels. -/
theorem exists_perm_coreEig (m : SpikedModelR μ n d rk) :
    ∃ σ : Equiv.Perm (Fin rk),
      (∀ k, m.toStack.coreEig (σ k) = m.θ k ^ 2) ∧
      ∀ (k : Fin rk) (N : ℕ),
        m.toStack.spikeVec (σ k) N = m.col N k
          ∨ m.toStack.spikeVec (σ k) N = -m.col N k := by
  refine ⟨corePerm rk, fun k => ?_, fun k N => ?_⟩
  · rw [m.coreEig_eq (corePerm rk k), corePerm, Equiv.apply_symm_apply]
  · have h := m.spikeVec_eq_or (corePerm rk k) N
    rwa [corePerm, Equiv.apply_symm_apply] at h

/-! ### 4. The sign cancels -/

/-- **Item C1.4.** The consequence Track C uses: `specProj` is linear and the statement squares
the norm, so the sign of the eigenvector never leaves this file. -/
theorem normSq_specProj_spikeVec_eq (m : SpikedModelR μ n d rk) {σ : Equiv.Perm (Fin rk)}
    (hσ : ∀ (k : Fin rk) (N : ℕ),
      m.toStack.spikeVec (σ k) N = m.col N k ∨ m.toStack.spikeVec (σ k) N = -m.col N k)
    (k : Fin rk) (N : ℕ) (ω : Ω N) (S : Set ℝ) :
    ‖specProj (m.toStack.gram N ω) S (m.toStack.spikeVec (σ k) N)‖ ^ 2
      = ‖specProj (m.toStack.gram N ω) S (m.col N k)‖ ^ 2 := by
  rcases hσ k N with h | h
  · rw [h]
  · rw [h, map_neg, norm_neg]

/-! ### 5. The tail shift -/

/-- **The shifted table.** The model reindexed by `N ↦ N + k`: same `θ`, every sequence
evaluated at `N + k`. The mirrors are `SpikedModel.shift` (`RMT/TailShift.lean:48`) and
`RankRStack.shift` (`RankR/RMT/ShiftR.lean:55`). -/
def shift (m : SpikedModelR μ n d rk) (k : ℕ) :
    SpikedModelR (fun N => μ (N + k)) (fun N => n (N + k)) (fun N => d (N + k)) rk where
  θ := m.θ
  U := fun N => m.U (N + k)
  V := fun N => m.V (N + k)
  Z := fun N => m.Z (N + k)
  hθnn := m.hθnn
  hθanti := m.hθanti
  hn := fun N => m.hn (N + k)
  hd := fun N => m.hd (N + k)
  hU := fun N => m.hU (N + k)
  hV := fun N => m.hV (N + k)
  hZ := fun N => m.hZ (N + k)

section ShiftBridges

variable (m : SpikedModelR μ n d rk) (k N : ℕ)

@[simp] theorem shift_θ : (m.shift k).θ = m.θ := rfl

@[simp] theorem shift_U : (m.shift k).U N = m.U (N + k) := rfl

@[simp] theorem shift_V : (m.shift k).V N = m.V (N + k) := rfl

@[simp] theorem shift_Z : (m.shift k).Z N = m.Z (N + k) := rfl

@[simp] theorem shift_E (ω : Ω (N + k)) : (m.shift k).E N ω = m.E (N + k) ω := rfl

@[simp] theorem shift_X (ω : Ω (N + k)) : (m.shift k).X N ω = m.X (N + k) ω := rfl

@[simp] theorem shift_col (l : Fin rk) : (m.shift k).col N l = m.col (N + k) l := rfl

@[simp] theorem shift_orthUnitR : (m.shift k).orthUnitR N = m.orthUnitR (N + k) := rfl

/-- The embedding commutes with the shift, so a `RankRStack` theorem proved for the shifted
stack applies to the shifted table with no transport. -/
theorem toStack_shift : (m.shift k).toStack = m.toStack.shift k := rfl

end ShiftBridges

/-- The regime is a tail property. -/
theorem shift_Regime {c : ℝ} (m : SpikedModelR μ n d rk) (k : ℕ) (hreg : m.Regime c) :
    (m.shift k).Regime c :=
  ⟨(tendsto_add_atTop_iff_nat (f := n) k).mpr hreg.1,
    (tendsto_add_atTop_iff_nat (f := d) k).mpr hreg.2.1,
    (tendsto_add_atTop_iff_nat (f := fun N => (n N : ℝ) / (d N : ℝ)) k).mpr hreg.2.2⟩

/-- The Gaussian noise law is stated index by index, so it shifts by evaluation. -/
theorem shift_GaussianNoise (m : SpikedModelR μ n d rk) (k : ℕ) (hG : m.GaussianNoise) :
    (m.shift k).GaussianNoise := fun N => hG (N + k)

/-- **The transfer.** A limit in probability of an index overlap for the shifted table gives
the limit for the table. Task C5.2 removes the side condition `∀ N, n N = rk + pp N` with
this lemma and `RankRStack.exists_shift_lt`. -/
theorem tendstoInProb_overlapIdx_of_shift (m : SpikedModelR μ n d rk) (k kk : ℕ)
    (l : Fin rk) {a : ℝ}
    (h : TendstoInProb (fun N => μ (N + k))
      (fun N ω => overlapIdx ((m.shift k).X N ω) kk ((m.shift k).col N l)) a) :
    TendstoInProb μ (fun N ω => overlapIdx (m.X N ω) kk (m.col N l)) a :=
  SpikedModel.tendstoInProb_of_shift k h

end SpikedModelR

end StackedSVD
