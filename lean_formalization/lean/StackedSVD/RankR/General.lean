/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.Defs
import StackedSVD.LinAlg.SpecIdx
import StackedSVD.LinAlg.SpecInvReindex
import StackedSVD.LinAlg.KyFan

/-!
# Section 7 at general `r_i`: model, block objects and statements

STATUS 2026-09-01: this file is `sorry` free. The five results of the old section 5 moved to
`RankR/GeneralMain.lean` on 2026-09-01 (task T6), because their proofs need `RankR/GramR.lean`
and `GramR.lean` imports this file. Every definition, the block bridge lemmas (`ABlock_intra`
and the four reductions of section 6), and `tableLawR_of_singleTableLaw` (task B5, proved
2026-09-01) stay here and are proved. `overlapIdx_zero` moved to `LinAlg/SpecIdx.lean` with
the index-spectral group and is proved there (task B1).
The review note is `notes/archive/rank_r_general.md` (task T8 of `notes/RANK_R_PLAN.md`).

`assum:unaligned` (`main_paper.tex:753`) gives table `i` an `r_i`-dimensional right singular
subspace `V R_i` inside one shared `r`-dimensional subspace `V`:

```
X_i = U_i Θ_i (V R_i)ᵀ + E_i,   R_i ∈ O(r, r_i),  Θ_i = diag(θ_i1, …, θ_i r_i),  U_i ∈ O(n_i, r_i).
```

`RankR/Defs.lean` takes the slice `r_i = 1`, where the per-table input is the proved rank-one
`SpikedModel.SingleTableLaw`. This file removes that restriction. Every object of
`RankR/Defs.lean` reappears with the table index `i` replaced by the block index `(i, j)`,
`i ∈ [M]`, `j ∈ [r_i]`, and `M` replaced by `r̃ = ∑_i r_i` (`main_paper.tex:765`).

## Content

1. `eigSetIdx`, `specProjIdx`, `overlapIdx`, `vEig`, `SimpleSpec`: the spectral projector at one
   **eigenvalue index**, not at the top eigenvalue. At index `0` they are `topProj`, `overlap`,
   `vMax` and `TopSimple` of `Defs.lean` and `SVDStack/Defs.lean`.
2. `SpikedModelR`: one table with `r_i` spikes. `TableLawR` is the rank-`r_i` twin of
   `SpikedModel.SingleTableLaw`, the only new hypothesis structure of this file.
3. `UnalignedModelR`: `M` such tables, a shared `V` with orthonormal columns and per-table
   `R_i` with orthonormal columns, with `(tbl i).V = V R_i`.
4. `rtot`, `blk`, `betaFlat`, `BBlock`, `DBlock`, `ABlock`: the `r̃ × r̃` block objects
   `B_R = diag(β) R_stack`, `D = diag(1 - β_ij²)` and `A_{β,R} = B_R B_Rᵀ + D`, with the paper's
   entrywise form `ABlock_apply`.
5. `perfRG`, `limitRG` and their weighted twins `perfRGW`, `limitRGW`, `optWG`, `limitOptG`.
6. The statements moved to `RankR/GeneralMain.lean`:
   `lem_general_rank_delocalization_general`, `gramR_general`, `VtV_tendsto_general`,
   `prop_general_rank_unweighted_svdstack_general`, `thm_gen_rank_weight_svdstak_general_r`.
7. The reductions to the `r_i = 1` slice: `tableLawR_of_singleTableLaw`,
   `BBlock_one_eq_BR`, `DBlock_one_eq_Dmat`, `ABlock_one_eq_AbetaR`, `limitRG_one_eq_limitR`.

## The index `(i, j)` and the index `p ∈ [r̃]`

The paper indexes the rows of `Ṽ`, of `B_R` and of `A_{β,R}` by the pair `(i, j)`. Lean needs
one `Fin` type, because `specInvTop`, `TopGap` and `Matrix.IsHermitian.eigenvalues₀` of
`LinAlg/SpecProjPerturb.lean` are stated for `Matrix (Fin p) (Fin p) ℝ`. The bridge is
`finSigmaFinEquiv : ((i : Fin M) × Fin (rk i)) ≃ Fin (∑ i, rk i)`, the same equivalence that
`MultiTableModel.stackZ` and `UnalignedModel.stackX` use for the stacked row index. `blk p` is
the pair of the flat index `p`, and every entry of every block object is read through it.

## The correspondence between the spike index and the eigenvalue index

`Ṽ` has one row per spike, `v̂_ij`, and `v̂_ij` must be the singular vector of `X_i` that tracks
the spike `θ_ij`. The model therefore orders the spikes of a table: `SpikedModelR.hθanti` says
`θ_i1 > θ_i2 > … > θ_i r_i`. This is the paper's own hypothesis at `main_paper.tex:768` ("the
diagonal entries of `Θ_i` are distinct"), made into an order, which is free: a simultaneous
permutation of the columns of `U_i`, `Θ_i` and `R_i` leaves the model unchanged. With it, the
`j`-th spike is the `j`-th largest eigenvalue of `X_iᵀ X_i` above the detection threshold, so
`v̂_ij := vEig (X_iᵀ X_i) j` is the right vector, and no permutation has to be carried.

`rk_le_d` and `norm_colVecG` moved here from `RankR/GeneralMain.lean`;
`tableLawR_of_singleTableLaw` and `SpikedModel.tableLawR_of_singleTableLaw'` moved out
to `RankR/GeneralMain.lean` (F28, 2026-09-08).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal MatrixOrder

namespace StackedSVD

/-! ### 0. The spectral projector at one eigenvalue index

`eigSetIdx`, `specProjIdx`, `overlapIdx`, `vEig`, `SimpleSpec` and `vEig_zero` moved to
`StackedSVD/LinAlg/SpecIdx.lean` on 2026-09-01 (task B1), together with the lemmas that use
them, `overlapIdx_zero` among them. This file consumes that interface. -/

/-! ### 1. The rank-`r_i` table and its law -/

/-- One table of `assum:unaligned` with `rk` spikes: `X = U Θ Vᵀ + d^{-1/2} Z` with `U` and `V`
of orthonormal columns and `Θ = diag(θ_1, …, θ_rk)`. The rank-`rk` twin of `SpikedModel`.

`hθanti` orders the spikes strictly, which is the paper's distinctness hypothesis
(`main_paper.tex:768`) plus the free choice of column order; see the header. `hθnn` is
`assum:unaligned`'s "diagonal with positive entries" weakened to `0 ≤ θ_k` (F8, 2026-09-05):
with `hθanti` it still forces at most one zero spike per table, at the last index
(`θ_eq_zero_unique`), and a supercritical spike is positive (`θ_pos_of_sup`). The
weakening matches `SpikedModel.hθ : 0 ≤ θ` at rank one and makes `SpikedModel.toRankR`
total. Neither field is read by any Layer 1 proof below; both are there because they are
what makes `TableLawR` satisfiable. -/
structure SpikedModelR {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    (μ : ∀ N, Measure (Ω N)) (n d : ℕ → ℕ) (rk : ℕ) where
  /-- the `rk` signal strengths, strictly decreasing -/
  θ : Fin rk → ℝ
  /-- left singular vectors, deterministic, orthonormal columns -/
  U : (N : ℕ) → Matrix (Fin (n N)) (Fin rk) ℝ
  /-- right singular vectors of this table, deterministic, orthonormal columns -/
  V : (N : ℕ) → Matrix (Fin (d N)) (Fin rk) ℝ
  /-- unscaled noise -/
  Z : (N : ℕ) → Ω N → Matrix (Fin (n N)) (Fin (d N)) ℝ
  hθnn : ∀ k, 0 ≤ θ k
  hθanti : StrictAnti θ
  hn : ∀ N, 0 < n N
  hd : ∀ N, 0 < d N
  hU : ∀ N, (U N)ᵀ * U N = 1
  hV : ∀ N, (V N)ᵀ * V N = 1
  hZ : ∀ N, Measurable (Z N)

namespace SpikedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n d : ℕ → ℕ} {rk : ℕ}

/-- Only the last spike may vanish: `θ_k > θ_{k+1} ≥ 0` for every non-final `k` (F8). -/
theorem θ_pos_of_ne_last (m : SpikedModelR μ n d rk) {k : Fin rk} (hk : (k : ℕ) + 1 < rk) :
    0 < m.θ k :=
  lt_of_le_of_lt (m.hθnn ⟨(k : ℕ) + 1, hk⟩) (m.hθanti (Fin.lt_def.mpr (Nat.lt_succ_self _)))

/-- A supercritical spike is positive. Mirror of the rank-one `theta_pos` of `RMT/R5.lean`. -/
theorem θ_pos_of_sup (m : SpikedModelR μ n d rk) {c : ℝ} (hc : 0 < c) {k : Fin rk}
    (hsup : c < m.θ k ^ 4) : 0 < m.θ k := by
  rcases (m.hθnn k).lt_or_eq with h | h
  · exact h
  · rw [← h] at hsup
    norm_num at hsup
    linarith

/-- At most one spike per table is zero. -/
theorem θ_eq_zero_unique (m : SpikedModelR μ n d rk) {k l : Fin rk} (hk : m.θ k = 0)
    (hl : m.θ l = 0) : k = l :=
  m.hθanti.injective (hk.trans hl.symm)

/-- Scaled noise `E_N = d_N^{-1/2} Z_N`. -/
noncomputable def E (m : SpikedModelR μ n d rk) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (n N)) (Fin (d N)) ℝ :=
  (Real.sqrt (d N))⁻¹ • m.Z N ω

/-- `X_N ω = U_N Θ V_Nᵀ + E_N ω`. -/
noncomputable def X (m : SpikedModelR μ n d rk) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (n N)) (Fin (d N)) ℝ :=
  m.U N * Matrix.diagonal m.θ * (m.V N)ᵀ + m.E N ω

/-- The `k`-th spike direction, the `k`-th column of `V_N`. -/
noncomputable def col (m : SpikedModelR μ n d rk) (N : ℕ) (k : Fin rk) :
    EuclideanSpace ℝ (Fin (d N)) :=
  WithLp.toLp 2 fun l => m.V N l k

/-- `assum:general_noise`, Gaussian case, as a law. -/
def GaussianNoise (m : SpikedModelR μ n d rk) : Prop :=
  ∀ N, HasLaw (m.Z N) (gaussianMatrix (n N) (d N)) (μ N)

/-- `eq:RMT_limit` for one table. -/
def Regime (_m : SpikedModelR μ n d rk) (c : ℝ) : Prop :=
  Tendsto n atTop atTop ∧ Tendsto d atTop atTop ∧
    Tendsto (fun N => (n N : ℝ) / d N) atTop (𝓝 c)

/-- Unit vectors orthogonal to the whole signal span of the table. The rank-`rk` twin of
`SpikedModel.orthUnit`; at `rk = 1` the two agree. -/
def orthUnitR (m : SpikedModelR μ n d rk) (N : ℕ) : Set (EuclideanSpace ℝ (Fin (d N))) :=
  {w | ‖w‖ = 1 ∧ ∀ k, ⟪w, m.col N k⟫_ℝ = 0}

/-- `prop:single_table` at rank `rk`, in projector form, one field per thing the proofs of
this file read. The rank-`rk` twin of `SpikedModel.SingleTableLaw` and the only new hypothesis
structure of Section 7 at general `r_i`.

* `align`: the `k`-th singular subspace of `X` overlaps the `k`-th spike direction by `β_k²`.
* `cross`: it does not overlap any **other** spike direction. This is where the distinctness of
  the `θ_k` (`SpikedModelR.hθanti`) is used; it is false at a tie.
* `delocUniform`: the `k`-th singular subspace does not overlap any direction orthogonal to the
  signal span, uniformly. This is the field the cross-table Fubini step consumes.
* `simple`: each of the top `rk` eigenvalues of `XᵀX` is simple almost surely, so each
  `specProjIdx` is a rank-one projector and every overlap is a squared inner product.

`SingleTableLaw.lamMax` has no twin here: `λ_k(XᵀX) → ρ²(θ_k, c)` is consumed by no statement
of this file (`CLAUDE.md` rule 5). It is a two-line addition when a consumer appears. -/
structure TableLawR (m : SpikedModelR μ n d rk) (c : ℝ) : Prop where
  align : ∀ k : Fin rk,
    TendstoInProb μ (fun N ω => overlapIdx (m.X N ω) (k : ℕ) (m.col N k)) (betaSq (m.θ k) c)
  cross : ∀ k l : Fin rk, k ≠ l →
    TendstoInProb μ (fun N ω => overlapIdx (m.X N ω) (k : ℕ) (m.col N l)) 0
  delocUniform : ∀ k : Fin rk, ∀ ε > 0,
    Tendsto (fun N => ⨆ w ∈ m.orthUnitR N, μ N {ω | ε ≤ overlapIdx (m.X N ω) (k : ℕ) w})
      atTop (𝓝 0)
  simple : ∀ N, ∀ᵐ ω ∂(μ N),
    SimpleSpec ((m.X N ω)ᵀ * m.X N ω) (isHermitian_transpose_mul_self (m.X N ω)) rk


/-- The model field `hV` forces `rk ≤ d N`: the `rk` columns of `V N` are independent in
`ℝ^{d N}`. This puts the eigenvalue index `j < rk` of a spike in range for `vEig` and
`specProjIdx`, both of which take a junk value out of range. -/
theorem rk_le_d (t : SpikedModelR μ n d rk) (N : ℕ) : rk ≤ d N := by
  have h1 : ((t.V N)ᵀ * t.V N).rank = rk := by
    rw [t.hV N]
    simp
  have h2 : ((t.V N)ᵀ * t.V N).rank ≤ (t.V N).rank := Matrix.rank_mul_le_right _ _
  have h3 : (t.V N).rank ≤ d N := by simpa using Matrix.rank_le_card_height (t.V N)
  omega

end SpikedModelR

/-! ### 2. The unaligned model at general `r_i` -/

/-- `assum:unaligned` (`main_paper.tex:753`) with no restriction on the `r_i`. `M` tables live
on one probability space per `N`. A shared `V N ∈ ℝ^{d_N × r}` has orthonormal columns, table
`i` carries a fixed `R_i ∈ O(r, r_i)`, and the right singular subspace of table `i` is
`V R_i`.

`Rank(∑_i R_i R_iᵀ) = r` of the paper is **not** a field, for the same reason as in
`UnalignedModel`: the unweighted proposition does not use it, and the weighted theorem takes
the stronger `rank B_R = r` as an explicit argument (decision D16 of `notes/FLAGGED.md`). -/
structure UnalignedModelR {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    (μ : ∀ N, Measure (Ω N)) (M : ℕ) (n : Fin M → ℕ → ℕ) (d : ℕ → ℕ) (r : ℕ)
    (rk : Fin M → ℕ) where
  /-- the `M` rank-`r_i` spiked tables -/
  tbl : (i : Fin M) → SpikedModelR μ (n i) d (rk i)
  /-- the shared ambient subspace, as a matrix with orthonormal columns -/
  V : (N : ℕ) → Matrix (Fin (d N)) (Fin r) ℝ
  /-- the alignment matrices `R_i ∈ O(r, r_i)`, fixed in `N` -/
  R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ
  hV : ∀ N, (V N)ᵀ * V N = 1
  hR : ∀ i, (R i)ᵀ * R i = 1
  hv : ∀ i N, (tbl i).V N = V N * R i

/-! ### 3. The block index and the block objects

`rtot rk` is the paper's `r̃ = ∑_i r_i` (`main_paper.tex:765`) and `blk p` is the pair `(i, j)`
of the flat index `p`. -/

/-- `r̃ = ∑_i r_i`, the number of spikes across all tables. -/
abbrev rtot {M : ℕ} (rk : Fin M → ℕ) : ℕ := ∑ i, rk i

/-- The block index `(i, j)` of a flat index `p ∈ [r̃]`. -/
def blk {M : ℕ} {rk : Fin M → ℕ} (p : Fin (rtot rk)) : (i : Fin M) × Fin (rk i) :=
  finSigmaFinEquiv.symm p

/-- The flat index of a block index `(i, j)`. -/
def flat {M : ℕ} {rk : Fin M → ℕ} (i : Fin M) (j : Fin (rk i)) : Fin (rtot rk) :=
  finSigmaFinEquiv ⟨i, j⟩

theorem blk_flat {M : ℕ} {rk : Fin M → ℕ} (i : Fin M) (j : Fin (rk i)) :
    blk (flat i j) = ⟨i, j⟩ := by
  rw [blk, flat, Equiv.symm_apply_apply]

section Block

variable {M r : ℕ} {rk : Fin M → ℕ}

/-- The vector `β ∈ ℝ^{r̃}` of the paper (`main_paper.tex:789`), read off the per-table family
`β_ij`. -/
noncomputable def betaFlat (β : (i : Fin M) → Fin (rk i) → ℝ) (p : Fin (rtot rk)) : ℝ :=
  β (blk p).1 (blk p).2

/-- `B_R = diag(β) R_stack ∈ ℝ^{r̃ × r}` (`main_paper.tex:1951`): row `(i, j)` is
`β_ij (R_i)_jᵀ`. -/
noncomputable def BBlock (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) : Matrix (Fin (rtot rk)) (Fin r) ℝ :=
  Matrix.of fun p k => betaFlat β p * R (blk p).1 k (blk p).2

/-- `D = diag(1 - β_ij²) ∈ ℝ^{r̃ × r̃}` (`main_paper.tex:2040`). -/
noncomputable def DBlock (β : (i : Fin M) → Fin (rk i) → ℝ) :
    Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ :=
  Matrix.diagonal fun p => 1 - betaFlat β p ^ 2

/-- `A_{β,R} = B_R B_Rᵀ + D` (`main_paper.tex:2049`). The paper's entrywise definition
(`main_paper.tex:779`) is `ABlock_apply`, `ABlock_diag` and `ABlock_intra`. -/
noncomputable def ABlock (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) :
    Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ :=
  BBlock β R * (BBlock β R)ᵀ + DBlock β

theorem isHermitian_ABlock (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) : (ABlock β R).IsHermitian :=
  Matrix.IsHermitian.add (isHermitian_mul_transpose_self (BBlock β R))
    (Matrix.isHermitian_diagonal _)

/-- The paper's entrywise definition of `A_{β,R}` (`main_paper.tex:779`):
`[A]_{(i,j),(i',j')} = β_ij β_{i'j'} (R_i)_jᵀ (R_{i'})_{j'}` off the diagonal, and `1` on it
once `ABlock_diag` is applied. -/
theorem ABlock_apply (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (p q : Fin (rtot rk)) :
    ABlock β R p q
      = betaFlat β p * betaFlat β q *
          ((R (blk p).1)ᵀ * R (blk q).1) (blk p).2 (blk q).2
        + (if p = q then 1 - betaFlat β p ^ 2 else 0) := by
  simp only [ABlock, Matrix.add_apply, Matrix.mul_apply, Matrix.transpose_apply, BBlock,
    Matrix.of_apply, DBlock, Matrix.diagonal_apply]
  congr 1
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun k _ => by ring

/-- The diagonal of `A_{β,R}` is `1`, which is the diagonal of `Ṽ Ṽᵀ` at every `N`. This is
where `R_i ∈ O(r, r_i)` is used. -/
theorem ABlock_diag {β : (i : Fin M) → Fin (rk i) → ℝ}
    {R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ} (hR : ∀ i, (R i)ᵀ * R i = 1)
    (p : Fin (rtot rk)) : ABlock β R p p = 1 := by
  rw [ABlock_apply, if_pos rfl, hR (blk p).1, Matrix.one_apply_eq]
  ring

/-- Inside one table the off-diagonal entries of `A_{β,R}` vanish, which is the paper's
`(Ṽ Ṽᵀ)_{(i,j),(i,j')} = 0` for `j ≠ j'` (`main_paper.tex:1972`). Again `R_i ∈ O(r, r_i)`. -/
theorem ABlock_intra {β : (i : Fin M) → Fin (rk i) → ℝ}
    {R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ} (hR : ∀ i, (R i)ᵀ * R i = 1)
    {p q : Fin (rtot rk)} (hpq : p ≠ q) (hsame : (blk p).1 = (blk q).1) :
    ABlock β R p q = 0 := by
  obtain ⟨i, j, rfl⟩ : ∃ i : Fin M, ∃ j : Fin (rk i), flat i j = p :=
    ⟨(blk p).1, (blk p).2, Equiv.apply_symm_apply finSigmaFinEquiv p⟩
  obtain ⟨i', j', rfl⟩ : ∃ i' : Fin M, ∃ j' : Fin (rk i'), flat i' j' = q :=
    ⟨(blk q).1, (blk q).2, Equiv.apply_symm_apply finSigmaFinEquiv q⟩
  simp only [blk_flat] at hsame
  subst hsame
  have hjne : j ≠ j' := fun h => hpq (congrArg (flat i) h)
  have hentry : ((R i)ᵀ * R i) j j' = 0 := by
    rw [hR i, Matrix.one_apply, if_neg hjne]
  rw [ABlock_apply, if_neg hpq]
  simp only [betaFlat, add_zero]
  rw [blk_flat, blk_flat, hentry]
  ring

/-- `A_{β,R} - D = B_R B_Rᵀ` is positive semidefinite, so `A_{β,R} ⪰ D`. -/
theorem ABlock_posSemidef_sub_DBlock (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) : (ABlock β R - DBlock β).PosSemidef := by
  have h : ABlock β R - DBlock β = BBlock β R * (BBlock β R)ᵀ := by
    rw [ABlock, add_sub_cancel_right]
  rw [h]
  simpa using Matrix.posSemidef_self_mul_conjTranspose (BBlock β R)

/-- `D ≻ 0` when every `β_ij` lies in `[0, 1)`. -/
theorem DBlock_posDef {β : (i : Fin M) → Fin (rk i) → ℝ} (h0 : ∀ i j, 0 ≤ β i j)
    (h1 : ∀ i j, β i j < 1) : (DBlock (rk := rk) β).PosDef := by
  refine Matrix.PosDef.diagonal fun p => ?_
  have ha := h0 (blk p).1 (blk p).2
  have hb := h1 (blk p).1 (blk p).2
  rw [betaFlat]
  nlinarith

/-- `A_{β,R} ≻ 0` when every `β_ij` lies in `[0, 1)`: the inverse inside
`specInvTop (A_{β,R}) r` never meets a zero eigenvalue. -/
theorem ABlock_posDef {β : (i : Fin M) → Fin (rk i) → ℝ}
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (h0 : ∀ i j, 0 ≤ β i j)
    (h1 : ∀ i j, β i j < 1) : (ABlock β R).PosDef := by
  rw [ABlock]
  exact posSemidef_add_posDef
    (by simpa using Matrix.posSemidef_self_mul_conjTranspose (BBlock β R))
    (DBlock_posDef h0 h1)

/-- The limit of `prop:general_rank_unweighted_svdstack` at general `r_i`, in the trace form of
`limitR`: `tr(B_Rᵀ (specInvTop A_{β,R} r) B_R)`, the paper's
`‖Λ^{-1/2} Qᵀ diag(β) R_stack‖_F²` (`main_paper.tex:802`). -/
noncomputable def limitRG (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) : ℝ :=
  Matrix.trace ((BBlock β R)ᵀ * specInvTop (ABlock β R) (isHermitian_ABlock β R) r *
    BBlock β R)

/-! #### The weighted block objects

`main_paper.tex:886` weights `Ṽ` by a full matrix `W ∈ ℝ^{r̃ × r̃}`, so `W` is not required to
be block diagonal. The optimum is the diagonal `W⋆ = D^{-1/2}` (`main_paper.tex:895`), whose
entries are the rank-one optimal weights `optW` of `SVDStack/Defs.lean`. -/

/-- `W A_{β,R} Wᵀ`, the limit of the weighted Gram matrix. -/
noncomputable def ABlockW (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) :
    Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ :=
  W * ABlock β R * Wᵀ

theorem isHermitian_ABlockW (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) : (ABlockW W β R).IsHermitian := by
  have h := Matrix.isHermitian_mul_mul_conjTranspose (A := ABlock β R) W (isHermitian_ABlock β R)
  rwa [Matrix.conjTranspose_eq_transpose_of_trivial] at h

/-- `W B_R`, the limit of `(W Ṽ) V`. -/
noncomputable def BBlockW (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) : Matrix (Fin (rtot rk)) (Fin r) ℝ :=
  W * BBlock β R

/-- The limit of the weighted performance at the weight matrix `W`. -/
noncomputable def limitRGW (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ)
    (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) : ℝ :=
  Matrix.trace ((BBlockW W β R)ᵀ *
    specInvTop (ABlockW W β R) (isHermitian_ABlockW W β R) r * BBlockW W β R)

/-- `W⋆ = D^{-1/2} = diag(1/√(1 - β_ij²))` (`main_paper.tex:895`). Its entries are the rank-one
optimal weights `optW`, which is the paper's statement that the optimal weighting does not
depend on the `R_i`. -/
noncomputable def optWG (β : (i : Fin M) → Fin (rk i) → ℝ) :
    Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ :=
  Matrix.diagonal (optW (betaFlat β))

/-- `A_{β,R}^{1/2}` from the continuous functional calculus. -/
noncomputable def ABlockSqrt (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) :
    Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ :=
  CFC.sqrt (ABlock β R)

theorem transpose_ABlockSqrt (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) :
    (ABlockSqrt β R)ᵀ = ABlockSqrt β R := by
  have h := (Matrix.nonneg_iff_posSemidef.mp (CFC.sqrt_nonneg (ABlock β R))).isHermitian
  rwa [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial] at h

theorem isHermitian_DBlock (β : (i : Fin M) → Fin (rk i) → ℝ) :
    (DBlock (rk := rk) β).IsHermitian :=
  Matrix.isHermitian_diagonal _

/-- `A_{β,R}^{-1/2} D A_{β,R}^{-1/2}`, the matrix whose bottom `r` eigenvalues give `L⋆`
(`main_paper.tex:898`). -/
noncomputable def DcongG (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) :
    Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ :=
  (ABlockSqrt β R)⁻¹ * DBlock β * (ABlockSqrt β R)⁻¹

theorem isHermitian_DcongG (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) : (DcongG β R).IsHermitian :=
  isHermitian_inv_mul_mul_inv (isHermitian_DBlock β) (transpose_ABlockSqrt β R)

/-- `L⋆ = r - ∑_{ℓ=1}^{r} λ_{r̃ + 1 - ℓ}(A_{β,R}^{-1/2} D A_{β,R}^{-1/2})`
(`main_paper.tex:898`). The paper's index `r̃ + 1 - ℓ` at `ℓ = k + 1` is `Fin.rev`, the index
form of `LinAlg/KyFan.lean`. -/
noncomputable def limitOptG (β : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ)
    (hr : r ≤ Fintype.card (Fin (rtot rk))) : ℝ :=
  (r : ℝ) - ∑ k : Fin r, (isHermitian_DcongG β R).eigenvalues₀ (Fin.castLE hr k).rev

end Block

/-! ### 4. The estimators on an `UnalignedModelR` -/

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- Joint law of all noise matrices at level `N`: independent Gaussian tables. The same
predicate as `MultiTableModel.JointGaussianNoise` and `UnalignedModel.JointGaussianNoise`. -/
def JointGaussianNoise (m : UnalignedModelR μ M n d r rk) : Prop :=
  ∀ N, HasLaw (fun ω i => (m.tbl i).Z N ω)
    (Measure.pi fun i : Fin M => gaussianMatrix (n i N) (d N)) (μ N)

/-- Gram matrix `X_iᵀ X_i` of table `i`. -/
noncomputable def tableGramG (m : UnalignedModelR μ M n d r rk) (i : Fin M) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (d N)) (Fin (d N)) ℝ := ((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω

theorem isHermitian_tableGramG (m : UnalignedModelR μ M n d r rk) (i : Fin M) (N : ℕ)
    (ω : Ω N) : (m.tableGramG i N ω).IsHermitian :=
  isHermitian_transpose_mul_self ((m.tbl i).X N ω)

/-- `v̂_ij`, the `j`-th right singular vector of table `i`, with the paper's sign convention
`⟪v̂_ij, (V R_i)_j⟫ ≥ 0` (`main_paper.tex:1958`). The spike order `hθanti` of the model is what
makes the sorted index `j` the right one; see the header. -/
noncomputable def vhatG (m : UnalignedModelR μ M n d r rk) (i : Fin M) (j : Fin (rk i))
    (N : ℕ) (ω : Ω N) : EuclideanSpace ℝ (Fin (d N)) :=
  if 0 ≤ ⟪vEig (m.tableGramG i N ω) (m.isHermitian_tableGramG i N ω) (j : ℕ),
      (m.tbl i).col N j⟫_ℝ then
    vEig (m.tableGramG i N ω) (m.isHermitian_tableGramG i N ω) (j : ℕ)
  else -vEig (m.tableGramG i N ω) (m.isHermitian_tableGramG i N ω) (j : ℕ)

/-- `Ṽ`, the `r̃ × d` matrix whose row `(i, j)` is `v̂_ijᵀ`. -/
noncomputable def VtG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (rtot rk)) (Fin (d N)) ℝ :=
  Matrix.of fun p l => m.vhatG (blk p).1 (blk p).2 N ω l

/-- `Ṽ Ṽᵀ`, the `r̃ × r̃` Gram matrix of the per-table estimates. -/
noncomputable def gramG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ := m.VtG N ω * (m.VtG N ω)ᵀ

theorem isHermitian_gramG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    (m.gramG N ω).IsHermitian :=
  isHermitian_mul_transpose_self (m.VtG N ω)

/-- `Ṽ V`, the `r̃ × r` matrix of the overlaps `⟪v̂_ij, V e_k⟫`. -/
noncomputable def VtVG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (rtot rk)) (Fin r) ℝ := m.VtG N ω * m.V N

/-- Column `k` of `V N`, as a vector of `EuclideanSpace`. -/
noncomputable def colVecG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (k : Fin r) :
    EuclideanSpace ℝ (Fin (d N)) := WithLp.toLp 2 fun l => m.V N l k


/-- The columns of `V N` are unit vectors. Mirror: `norm_colVec`. -/
theorem norm_colVecG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (k : Fin r) :
    ‖m.colVecG N k‖ = 1 := by
  have hVV := m.hV N
  have hsq : ‖m.colVecG N k‖ ^ 2 = 1 := by
    rw [← real_inner_self_eq_norm_sq, real_inner_eq_dotProduct]
    change ∑ l, (m.V N) l k * (m.V N) l k = 1
    have h : ∑ l, (m.V N) l k * (m.V N) l k = ((m.V N)ᵀ * m.V N) k k := by
      rw [Matrix.mul_apply]
      exact Finset.sum_congr rfl fun l _ => by rw [Matrix.transpose_apply]
    rw [h, hVV]
    simp
  nlinarith [norm_nonneg (m.colVecG N k)]

/-- Performance of svdstack at general `r_i`: `tr((Ṽ V)ᵀ (specInvTop (Ṽ Ṽᵀ) r) (Ṽ V))`, which
equals `‖V̂_svdstackᵀ V‖_F²` at any top-`r` eigenframe of `Ṽ Ṽᵀ` with nonnegative eigenvalues
(`main_paper.tex:1982`; `frobSq_vhatSvdstackG` in `RankR/GeneralFrob.lean`, which needs the
eigengap `λ_r > λ_{r+1}` for the frame to be the top-`r` one, and `λ_r > 0` for the inverse;
Track B audit C4). The identity was measured to 3.1e-15 over 96 draws at `M = 2`,
`r_i = (2, 1)` (`notes/archive/rank_r_general.md`, scan row M2). -/
noncomputable def perfRG (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) : ℝ :=
  Matrix.trace ((m.VtVG N ω)ᵀ *
    specInvTop (m.gramG N ω) (m.isHermitian_gramG N ω) r * m.VtVG N ω)

/-- `Ṽ_W = W Ṽ`, whose top `r` right singular vectors are `V̂_svdstack(W)`
(`main_paper.tex:886`). -/
noncomputable def VtWG (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (rtot rk)) (Fin (d N)) ℝ := W * m.VtG N ω

noncomputable def gramWG (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ := m.VtWG W N ω * (m.VtWG W N ω)ᵀ

theorem isHermitian_gramWG (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (N : ℕ) (ω : Ω N) :
    (m.gramWG W N ω).IsHermitian :=
  isHermitian_mul_transpose_self (m.VtWG W N ω)

noncomputable def VtVWG (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (rtot rk)) (Fin r) ℝ := m.VtWG W N ω * m.V N

/-- Performance of weighted svdstack at general `r_i`. At `W = 1` it is `perfRG`. -/
noncomputable def perfRGW (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (N : ℕ) (ω : Ω N) : ℝ :=
  Matrix.trace ((m.VtVWG W N ω)ᵀ *
    specInvTop (m.gramWG W N ω) (m.isHermitian_gramWG W N ω) r * m.VtVWG W N ω)

end UnalignedModelR

/-! ### 6. The reduction to the `r_i = 1` slice

At `rk = fun _ => 1` the block index `(i, j)` is the table index `i`, `r̃ = M`, and the objects
of this file are the objects of `RankR/Defs.lean`. The flat index type `Fin (∑ i : Fin M, 1)`
is not definitionally `Fin M`, so the statements read the block objects along the injection
`i ↦ flat i 0`, which is a bijection here. -/

section Reduce

variable {M r : ℕ}

/-- The alignment matrix of a unit vector, the `r_i = 1` shape of `R_i ∈ O(r, 1)`. -/
noncomputable def Rone (R : Fin M → EuclideanSpace ℝ (Fin r)) (i : Fin M) :
    Matrix (Fin r) (Fin ((fun _ : Fin M => 1) i)) ℝ :=
  Matrix.of fun k _ => R i k

/-- `i ↦ flat i 0` is injective at `rk = fun _ => 1`: both indices land in `Fin 1`, so the flat
index remembers only the table. Used by the three reductions below to move the diagonal
`if`-condition from the flat index back to the table index. -/
private theorem flat_zero_eq_iff {M : ℕ} {i i' : Fin M} :
    flat (rk := fun _ : Fin M => 1) i 0 = flat i' 0 ↔ i = i' := by
  constructor
  · intro h
    have h2 := congrArg (blk (rk := fun _ : Fin M => 1)) h
    rw [blk_flat, blk_flat] at h2
    exact congrArg Sigma.fst h2
  · rintro rfl; rfl

/-- `B_R` at `r_i = 1` is `BR` of `RankR/Defs.lean`. -/
theorem BBlock_one_eq_BR (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    (BBlock (rk := fun _ => 1) (fun i _ => β i) (Rone R)).submatrix
        (fun i : Fin M => flat i 0) id = BR β R := by
  ext i k
  simp [Matrix.submatrix_apply, BBlock, betaFlat, blk_flat, Rone, BR, Matrix.of_apply]

/-- `D` at `r_i = 1` is `Dmat` of `RankR/Defs.lean`. -/
theorem DBlock_one_eq_Dmat (β : Fin M → ℝ) :
    (DBlock (rk := fun _ : Fin M => 1) (fun i _ => β i)).submatrix
        (fun i : Fin M => flat i 0) (fun i : Fin M => flat i 0) = Dmat β := by
  ext i i'
  simp [Matrix.submatrix_apply, DBlock, Matrix.diagonal_apply, Dmat, betaFlat, blk_flat,
    flat_zero_eq_iff]

/-- `A_{β,R}` at `r_i = 1` is `AbetaR` of `RankR/Defs.lean` (`main_paper.tex:810`). -/
theorem ABlock_one_eq_AbetaR (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    (ABlock (rk := fun _ => 1) (fun i _ => β i) (Rone R)).submatrix
        (fun i : Fin M => flat i 0) (fun i : Fin M => flat i 0) = AbetaR β R := by
  ext i i'
  have hinner : ((Rone R i)ᵀ * Rone R i') 0 0 = ⟪R i, R i'⟫_ℝ := by
    rw [real_inner_eq_dotProduct, Matrix.mul_apply]
    simp [Rone, dotProduct, Matrix.transpose_apply]
  simp only [Matrix.submatrix_apply]
  rw [ABlock_apply, abetaR_apply]
  simp only [betaFlat]
  rw [blk_flat, blk_flat, hinner]
  simp [flat_zero_eq_iff]

/-- The limit at `r_i = 1` is `limitR` of `RankR/Defs.lean`. Unlike the three matrix
reductions above, this one needs the invariance of `tr(Bᵀ (specInvTop A r) B)` under a
simultaneous reindexing of `A` and `B` by a bijection of the index type. That invariance is
`trace_specInvTop_congr_submatrix` of `LinAlg/SpecInvReindex.lean`, and the bijection here is
`i ↦ flat i 0`, which is injective by `blk_flat` and surjective by a count of the two index
types. -/
theorem limitRG_one_eq_limitR (β : Fin M → ℝ) (R : Fin M → EuclideanSpace ℝ (Fin r)) :
    limitRG (rk := fun _ => 1) (fun i _ => β i) (Rone R) = limitR β R := by
  obtain ⟨e, he⟩ : ∃ e : Fin M ≃ Fin (rtot fun _ : Fin M => 1),
      ∀ i, e i = flat (rk := fun _ : Fin M => 1) i 0 := by
    refine ⟨Equiv.ofBijective (fun i : Fin M => flat (rk := fun _ : Fin M => 1) i 0) ?_,
      fun i => rfl⟩
    refine (Fintype.bijective_iff_injective_and_card _).mpr ⟨?_, by simp [rtot]⟩
    intro i i' h
    have h2 : blk (rk := fun _ : Fin M => 1) (flat i 0) = blk (flat i' 0) := congrArg blk h
    rw [blk_flat, blk_flat] at h2
    exact congrArg Sigma.fst h2
  have hcoe : (⇑e) = fun i : Fin M => flat (rk := fun _ : Fin M => 1) i 0 := funext he
  simp only [limitRG, limitR]
  refine (trace_specInvTop_congr_submatrix _ _ _ e _ _ _ ?_ ?_ r).symm
  · rw [hcoe]
    exact (ABlock_one_eq_AbetaR β R).symm
  · rw [hcoe]
    exact (BBlock_one_eq_BR β R).symm

end Reduce

/-! ### 7. `SpikedModel` as a one-spike `SpikedModelR`

`SpikedModel` (`Defs.lean`) is the rank-one model. `SpikedModelR μ n d 1` above is its
rank-`rk` twin at `rk = 1`. Both ask only `0 ≤ θ` (F8, 2026-09-05), so the map from
`SpikedModel` to `SpikedModelR` is total: `SpikedModel.toRankR`. The reverse map,
`SpikedModelR.toSpiked`, is in `RankR/GeneralFrob.lean`. With `toRankR`,
`tableLawR_of_singleTableLaw` above discharges `TableLawR c` for the constructed model
straight from `SingleTableLaw c`, in `tableLawR_of_singleTableLaw'`. -/

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {n d : ℕ → ℕ}

/-- A rank-one `SpikedModel` as a one-spike `SpikedModelR` (F8, 2026-09-05): total, because
both models ask only `0 ≤ θ`. -/
noncomputable def toRankR (m : SpikedModel μ n d) : SpikedModelR μ n d 1 where
  θ := fun _ => m.θ
  U N := Matrix.of fun i _ => m.u N i
  V N := Matrix.of fun l _ => m.v N l
  Z := m.Z
  hθnn := fun _ => m.hθ
  hθanti := Subsingleton.strictAnti _
  hn := m.hn
  hd := m.hd
  hU N := by
    have h : ∑ i, m.u N i ^ 2 = 1 := by
      have hs := EuclideanSpace.real_norm_sq_eq (m.u N)
      rw [m.hu N, one_pow] at hs
      exact hs.symm
    ext i j
    fin_cases i
    fin_cases j
    simp only [Matrix.mul_apply, Matrix.transpose_apply, Matrix.of_apply, Matrix.one_apply_eq]
    simpa [pow_two] using h
  hV N := by
    have h : ∑ l, m.v N l ^ 2 = 1 := by
      have hs := EuclideanSpace.real_norm_sq_eq (m.v N)
      rw [m.hv N, one_pow] at hs
      exact hs.symm
    ext i j
    fin_cases i
    fin_cases j
    simp only [Matrix.mul_apply, Matrix.transpose_apply, Matrix.of_apply, Matrix.one_apply_eq]
    simpa [pow_two] using h
  hZ := m.hZ

/-- The two models carry the same data matrix. -/
theorem toRankR_X (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) :
    m.toRankR.X N ω = m.X N ω := by
  ext i l
  simp [SpikedModelR.X, SpikedModelR.E, SpikedModel.X, SpikedModel.E, toRankR,
    Matrix.mul_apply, Matrix.vecMulVec_apply, mul_comm, mul_assoc]

/-- The one column of `V N` in `toRankR` is the spike direction `v N`. -/
theorem toRankR_col (m : SpikedModel μ n d) (N : ℕ) : m.toRankR.col N 0 = m.v N := rfl

/-- The one spike of `toRankR` is `m.θ`. -/
theorem toRankR_θ (m : SpikedModel μ n d) : m.toRankR.θ 0 = m.θ := rfl

/-- The regime is a statement about `n` and `d` only, so it transfers by `rfl`, as for
`SpikedModelR.toSpiked_Regime` (`RankR/GeneralFrob.lean`). -/
theorem toRankR_Regime (m : SpikedModel μ n d) (c : ℝ) : m.toRankR.Regime c ↔ m.Regime c :=
  Iff.rfl

/-- The two models have the same noise function, so Gaussian noise transfers by `rfl`, as for
`SpikedModelR.toSpiked_GaussianNoise`. -/
theorem toRankR_GaussianNoise (m : SpikedModel μ n d) :
    m.toRankR.GaussianNoise ↔ m.GaussianNoise :=
  Iff.rfl

end SpikedModel

end StackedSVD
