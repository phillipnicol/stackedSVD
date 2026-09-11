/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.StackSVDWeighted
import StackedSVD.RankR.SubspaceG
import StackedSVD.RankR.Aligned
import StackedSVD.LinAlg.SpecIdxPerturb
import StackedSVD.RankR.Subspace

/-!
# Rank-`r` weighted stacksvd: the scalars and the per-component stack

Track D item D1 of `notes/RANK_R_PLAN.md`. The review note is
`notes/archive/rankr_D1_statement.md`; this file carries the definitions and the facts that need no
probability limit. The limit theorem `thm:rank_r_stacksvd` (`main_paper.tex:2337`) is item D2
and is not stated here.

## The paper

In the exactly aligned rank-`r` model `X_i = ∑_j θ_ij u_ij v_jᵀ + E_i`
(`eq:rank_r_model`, `main_paper.tex:2298`) the paper builds one weighted stack per component
`j` (`eq:stacksvd_app_wij` to `eq:stacksvd_gammak`):

```
w_ij = θ_ij / √(θ_ij² + c_i),        X_stack^(j) = [w_1j X_1; …; w_Mj X_M],
θ̃_jk² = ∑_i θ_ij² θ_ik² / (θ_ij² + c_i),
γ_j   = the root x ∈ (0,1) of ∑_i θ_ij⁴ (1 - x)/(c_i + x θ_ij²) = 1.
```

The estimator `v̂_{j,stacksvd}` is the `ℓ_j`-th right singular vector of `X_stack^(j)`, where
`ℓ_j` is the rank of `θ̃_jj` among `{θ̃_jk}_k` in decreasing order (line
`alg_line:stacksvd_rank` of `alg:rank_r_stacksvd`).

## The four definitions here

* `Scalars.wStackR θ c j = optWstack (θ · j) c`, the paper's `w_·j`.
* `Scalars.thetaTildeSq θ c j k`, the paper's `θ̃_jk²`.
* `Scalars.gammaR θ c j = stackSVDLimitW (θ · j) c`, the paper's `γ_j`.
* `Scalars.ellR θ c j`, the paper's `ℓ_j`, as a **0-based** index: the number of components
  whose strength under the `j`-th weighting is strictly above `θ̃_jj`.

`γ_j` is **total**: `Scalars.stackSVDLimitW` returns `0` when no root exists, that is when
`∑_i θ_ij⁴/c_i ≤ 1`. This is `paper_edits.md` E3, and it matches the rank-1 weighted result
`thm:stacksvd_weighted`, whose limit is also `0` below the same threshold
(`Scalars.stackSVDLimitW_eq_zero`).

## The model side

`UnalignedModelR.stackXW m w` is `stackXG` with block row `i` scaled by `w i`, for any `rk`;
at `w = 1` it is `stackXG` (`stackXW_one`). On the exactly aligned family
`rk = alignedRk M r` it specializes to `stackXJ c j = X_stack^(j)`, with Gram `stackGramJ`,
selected eigenvector `vhatStackR` (sign fixed as `vhatG` does), squared overlap
`stackOverlapJ` and Frobenius aggregate `frobSqStackR`.

`stackOverlapJ c j N ω = overlapIdx (X_stack^(j)) ℓ_j v_j` is the projector form of the
paper's `(v_jᵀ v̂_j)²`; `stackOverlapJ_eq_inner_sq` is the identity with the paper's inner
product, under `SimpleIdx` at the `ℓ_j`-th eigenvalue (one index, not the whole top block).

No `sorry`, no new axiom.

`structure HeteroLawR` moved here from `RankR/StackMain.lean` (F28, 2026-09-08).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

/-- An orthogonal projector is a contraction, so a squared overlap at any sorted index is at
most `‖w‖²`. Mirror: `overlap_le_norm_sq` (`Spectral.lean`), which states the same bound at
index `0`. Private: it belongs in `LinAlg/SpecIdx.lean`, and only this file needs it today. -/
private theorem overlapIdx_le_norm_sq {p q : ℕ} (X : Matrix (Fin q) (Fin p) ℝ) (k : ℕ)
    (w : EuclideanSpace ℝ (Fin p)) : overlapIdx X k w ≤ ‖w‖ ^ 2 := by
  have h : ‖specProjIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) k w‖ ≤ ‖w‖ :=
    Submodule.norm_starProjection_apply_le _ w
  have h0 : (0 : ℝ) ≤ ‖specProjIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) k w‖ :=
    norm_nonneg _
  calc overlapIdx X k w
      = ‖specProjIdx (Xᵀ * X) (isHermitian_transpose_mul_self X) k w‖ ^ 2 := rfl
    _ ≤ ‖w‖ ^ 2 := by nlinarith

/-! ### 1. The scalars of `thm:rank_r_stacksvd` -/

namespace Scalars

variable {M r : ℕ}

/-- `w_ij = θ_ij/√(θ_ij² + c_i)` (`eq:stacksvd_app_wij`): the rank-1 optimal weights
`optWstack` applied to column `j` of the strength table. -/
noncomputable def wStackR (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (j : Fin r) : Fin M → ℝ :=
  optWstack (fun i => θ i j) c

theorem wStackR_apply (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (j : Fin r) (i : Fin M) :
    wStackR θ c j i = θ i j / Real.sqrt (θ i j ^ 2 + c i) := rfl

/-- `w_ij² = θ_ij²/(θ_ij² + c_i)`. -/
theorem wStackR_sq {θ : Fin M → Fin r → ℝ} {c : Fin M → ℝ} (hc : ∀ i, 0 < c i) (j : Fin r)
    (i : Fin M) : wStackR θ c j i ^ 2 = θ i j ^ 2 / (θ i j ^ 2 + c i) := by
  have hpos : 0 < θ i j ^ 2 + c i := add_pos_of_nonneg_of_pos (sq_nonneg _) (hc i)
  rw [wStackR_apply, div_pow, Real.sq_sqrt hpos.le]

/-- `θ̃_jk² = ∑_i θ_ij² θ_ik²/(θ_ij² + c_i)` (`eq:stacksvd_apptildTheta`): the strength of
component `k` under the `j`-th weighting. -/
noncomputable def thetaTildeSq (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (j k : Fin r) : ℝ :=
  ∑ i, θ i j ^ 2 * θ i k ^ 2 / (θ i j ^ 2 + c i)

/-- The paper's reading of `θ̃_jk²` (`main_paper.tex:2376`): `∑_i w_ij² θ_ik²`. -/
theorem thetaTildeSq_eq_sum_wStackR_sq {θ : Fin M → Fin r → ℝ} {c : Fin M → ℝ}
    (hc : ∀ i, 0 < c i) (j k : Fin r) :
    thetaTildeSq θ c j k = ∑ i, wStackR θ c j i ^ 2 * θ i k ^ 2 := by
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [wStackR_sq hc j i, div_mul_eq_mul_div]

/-- The diagonal strength `θ̃_jj² = ∑_i θ_ij⁴/(θ_ij² + c_i)`. -/
theorem thetaTildeSq_diag (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (j : Fin r) :
    thetaTildeSq θ c j j = ∑ i, θ i j ^ 4 / (θ i j ^ 2 + c i) := by
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [show θ i j ^ 2 * θ i j ^ 2 = θ i j ^ 4 by ring]

theorem thetaTildeSq_nonneg {θ : Fin M → Fin r → ℝ} {c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (j k : Fin r) : 0 ≤ thetaTildeSq θ c j k :=
  Finset.sum_nonneg fun i _ =>
    div_nonneg (by positivity) (add_pos_of_nonneg_of_pos (sq_nonneg _) (hc i)).le

/-- `γ_j` (`eq:stacksvd_gammak`), **total**: the unique root in `(0,1)` of
`∑_i θ_ij⁴(1-x)/(c_i + xθ_ij²) = 1` when `∑_i θ_ij⁴/c_i > 1`, and `0` otherwise
(`paper_edits.md` E3). This is the rank-1 limit `stackSVDLimitW` read on column `j`. -/
noncomputable def gammaR (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (j : Fin r) : ℝ :=
  stackSVDLimitW (fun i => θ i j) c

theorem gammaR_nonneg (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (j : Fin r) :
    0 ≤ gammaR θ c j :=
  stackSVDLimitW_nonneg _ _

/-- Below the detectability threshold `γ_j = 0`. -/
theorem gammaR_eq_zero {θ : Fin M → Fin r → ℝ} {c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    {j : Fin r} (h : ∑ i, θ i j ^ 4 / c i ≤ 1) : gammaR θ c j = 0 :=
  stackSVDLimitW_eq_zero hc h

/-- Above the threshold `γ_j` is the root of the paper's secular equation, and it lies in
`(0,1)`. -/
theorem gammaR_spec {θ : Fin M → Fin r → ℝ} {c : Fin M → ℝ} (hc : ∀ i, 0 < c i) {j : Fin r}
    (hthr : 1 < ∑ i, θ i j ^ 4 / c i) :
    gammaR θ c j ∈ Set.Ioo (0 : ℝ) 1 ∧ gW (fun i => θ i j) c (gammaR θ c j) = 1 :=
  stackSVDLimitW_spec (existsUnique_root hc hthr)

/-- `γ_j < 1` in both regimes. -/
theorem gammaR_lt_one {θ : Fin M → Fin r → ℝ} {c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (j : Fin r) : gammaR θ c j < 1 := by
  by_cases h : 1 < ∑ i, θ i j ^ 4 / c i
  · exact (gammaR_spec hc h).1.2
  · rw [gammaR_eq_zero hc (not_lt.mp h)]
    norm_num

/-- Above the threshold `γ_j > 0`. -/
theorem gammaR_pos {θ : Fin M → Fin r → ℝ} {c : Fin M → ℝ} (hc : ∀ i, 0 < c i) {j : Fin r}
    (hthr : 1 < ∑ i, θ i j ^ 4 / c i) : 0 < gammaR θ c j :=
  (gammaR_spec hc hthr).1.1

/-- `ℓ_j` of `alg:rank_r_stacksvd`, as a **0-based** sorted index: the number of components
`k` whose strength `θ̃_jk²` under the `j`-th weighting is strictly above `θ̃_jj²`. The paper's
1-based rank is `ellR + 1`. `eigenvalues₀` is antitone, so `ellR` is the index that
`overlapIdx` consumes. -/
noncomputable def ellR (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (j : Fin r) : ℕ :=
  (Finset.univ.filter fun k => thetaTildeSq θ c j j < thetaTildeSq θ c j k).card

/-- `ℓ_j < r`: the index `j` itself never enters the filter. -/
theorem ellR_lt (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (j : Fin r) : ellR θ c j < r := by
  have hj : j ∉ Finset.univ.filter fun k => thetaTildeSq θ c j j < thetaTildeSq θ c j k := by
    simp
  have hss : (Finset.univ.filter fun k => thetaTildeSq θ c j j < thetaTildeSq θ c j k) ⊂
      (Finset.univ : Finset (Fin r)) :=
    ⟨Finset.filter_subset _ _, fun h => hj (h (Finset.mem_univ j))⟩
  have := Finset.card_lt_card hss
  simpa [ellR] using this

/-- When `θ̃_jj²` is the strict maximum, the estimator is the top singular vector: `ℓ_j = 0`.
This is the paper's "if `θ_ij` follow the same ordering in each table then `ℓ_j = j`" at
`j = 0`. -/
theorem ellR_eq_zero_of_max {θ : Fin M → Fin r → ℝ} {c : Fin M → ℝ} {j : Fin r}
    (h : ∀ k, k ≠ j → thetaTildeSq θ c j k < thetaTildeSq θ c j j) : ellR θ c j = 0 := by
  rw [ellR, Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  intro k _
  rcases eq_or_ne k j with rfl | hk
  · exact lt_irrefl _
  · exact not_lt.mpr (h k hk).le

/-- When every table orders its spikes the same way and one table has a positive strength at
the weighting component `j`, the strength under weighting `j` is strictly decreasing in `k`.
This is the paper's sentence "if `θ_ij` follow the same ordering in each table, then the
ordering is preserved for each weighting" (`main_paper.tex:2370`). -/
theorem thetaTildeSq_strictAnti {θ : Fin M → Fin r → ℝ} {c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hnn : ∀ i k, 0 ≤ θ i k) (hanti : ∀ i, StrictAnti (θ i)) {j : Fin r}
    (hex : ∃ i, 0 < θ i j) : StrictAnti (thetaTildeSq θ c j) := by
  intro k l hkl
  obtain ⟨i₀, hi₀⟩ := hex
  refine Finset.sum_lt_sum (fun i _ => ?_) ⟨i₀, Finset.mem_univ _, ?_⟩
  · have hden : 0 < θ i j ^ 2 + c i := add_pos_of_nonneg_of_pos (sq_nonneg _) (hc i)
    have hlt : θ i l < θ i k := hanti i hkl
    have hsq : θ i l ^ 2 ≤ θ i k ^ 2 := by nlinarith [hnn i l]
    exact div_le_div_of_nonneg_right (mul_le_mul_of_nonneg_left hsq (sq_nonneg _)) hden.le
  · have hden : 0 < θ i₀ j ^ 2 + c i₀ := add_pos_of_nonneg_of_pos (sq_nonneg _) (hc i₀)
    have hlt : θ i₀ l < θ i₀ k := hanti i₀ hkl
    have hsq : θ i₀ l ^ 2 < θ i₀ k ^ 2 := by nlinarith [hnn i₀ l]
    have hj2 : (0 : ℝ) < θ i₀ j ^ 2 := pow_pos hi₀ 2
    rw [div_lt_div_iff₀ hden hden]
    nlinarith [mul_pos (mul_pos hj2 (sub_pos.mpr hsq)) hden]

/-- **`ℓ_j = j` whenever the spike order is the same in every table.** The paper's remark at
`main_paper.tex:2370`, and the reason the index machinery is only needed for unordered
`θ_ij`. -/
theorem ellR_eq_val {θ : Fin M → Fin r → ℝ} {c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hnn : ∀ i k, 0 ≤ θ i k) (hanti : ∀ i, StrictAnti (θ i)) {j : Fin r}
    (hex : ∃ i, 0 < θ i j) : ellR θ c j = (j : ℕ) := by
  have hsa : StrictAnti (thetaTildeSq θ c j) := thetaTildeSq_strictAnti hc hnn hanti hex
  have hset : (Finset.univ.filter fun k => thetaTildeSq θ c j j < thetaTildeSq θ c j k)
      = Finset.Iio j := by
    ext k
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_Iio]
    exact StrictAnti.lt_iff_gt hsa
  rw [ellR, hset, Fin.card_Iio]

/-! #### The paper's worked example (`main_paper.tex:2374`)

`θ_1 = (2, 1)`, `θ_2 = (1, 10)`, `c = (1, 1)`. Under the first weighting the strength of `v_1`
is `16/5 + 1/2 = 3.7` and the strength of `v_2` is `4/5 + 100/2 = 50.8`, so `ℓ_1 = 1` in the
0-based index (the paper's second position). -/

example : thetaTildeSq (M := 2) (r := 2) ![![2, 1], ![1, 10]] ![1, 1] 0 0 = 3.7 := by
  norm_num [thetaTildeSq, Fin.sum_univ_two]

example : thetaTildeSq (M := 2) (r := 2) ![![2, 1], ![1, 10]] ![1, 1] 0 1 = 50.8 := by
  norm_num [thetaTildeSq, Fin.sum_univ_two]

example : ellR (M := 2) (r := 2) ![![2, 1], ![1, 10]] ![1, 1] 0 = 1 := by
  have h0 : thetaTildeSq (M := 2) (r := 2) ![![2, 1], ![1, 10]] ![1, 1] 0 0 = 3.7 := by
    norm_num [thetaTildeSq, Fin.sum_univ_two]
  have h1 : thetaTildeSq (M := 2) (r := 2) ![![2, 1], ![1, 10]] ![1, 1] 0 1 = 50.8 := by
    norm_num [thetaTildeSq, Fin.sum_univ_two]
  rw [ellR, Finset.card_filter, Fin.sum_univ_two]
  simp only [h0, h1]
  norm_num

end Scalars

/-! ### 2. The weighted stack at general `r_i` -/

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- The weighted stack `[w_1 X_1; …; w_M X_M]` (`eq:stacksvd_appXstack`), with rows reindexed
from the block index `(i : Fin M) × Fin (n i N)` exactly as `stackXG` does. The rank-1 mirror
is `MultiTableModel.stackW`. -/
noncomputable def stackXW (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
  Matrix.reindex finSigmaFinEquiv (Equiv.refl (Fin (d N)))
    (Matrix.of fun (p : (i : Fin M) × Fin (n i N)) k => w p.1 * (m.tbl p.1).X N ω p.2 k)

/-- Block `i` of the weighted stack is `w_i X_i`. -/
theorem stackXW_apply (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N)
    (i : Fin M) (a : Fin (n i N)) (k : Fin (d N)) :
    m.stackXW w N ω (finSigmaFinEquiv ⟨i, a⟩) k = w i * (m.tbl i).X N ω a k := by
  simp [stackXW, Matrix.reindex_apply, Matrix.submatrix_apply]

/-- The weighted stack read at an arbitrary row index. -/
theorem stackXW_apply' (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N)
    (q : Fin (∑ i, n i N)) (k : Fin (d N)) :
    m.stackXW w N ω q k
      = w (finSigmaFinEquiv.symm q).1 *
        (m.tbl (finSigmaFinEquiv.symm q).1).X N ω (finSigmaFinEquiv.symm q).2 k := rfl

/-- At unit weights the weighted stack is `stackXG`. -/
theorem stackXW_one (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    m.stackXW 1 N ω = m.stackXG N ω := by
  ext q k
  rw [stackXW_apply', stackXG_apply']
  exact one_mul _

/-- Gram matrix of the weighted stack. -/
noncomputable def stackGramW (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : Matrix (Fin (d N)) (Fin (d N)) ℝ :=
  (m.stackXW w N ω)ᵀ * m.stackXW w N ω

theorem isHermitian_stackGramW (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : (m.stackGramW w N ω).IsHermitian :=
  isHermitian_transpose_mul_self _

/-- `X_stack(w)ᵀ X_stack(w) = ∑_i w_i² X_iᵀ X_i`. The weighted twin of `stackGramG_eq_sum`. -/
theorem stackGramW_eq_sum (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) :
    m.stackGramW w N ω
      = ∑ i, w i ^ 2 • (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω) := by
  ext k l
  simp only [stackGramW, Matrix.mul_apply, Matrix.transpose_apply, Matrix.sum_apply,
    Matrix.smul_apply, smul_eq_mul, Finset.mul_sum, stackXW_apply']
  rw [sum_stack_index (ν := fun i => n i N)
    (fun i a => w i * (m.tbl i).X N ω a k * (w i * (m.tbl i).X N ω a l))]
  exact Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun a _ => by ring

/-! ### 3. The exactly aligned family: `X_stack^(j)` and `(v_jᵀ v̂_j)²`

`rk = alignedRk M r` and `R_i = 1`, so spike `j` of table `i` is the shared `v_j` and the
strength table is `θ_ij`. -/

section AlignedStack

/-- The strength table `θ_ij` of the exactly aligned rank-`r` model. -/
noncomputable def thetaAligned (m : UnalignedModelR μ M n d r (alignedRk M r)) :
    Fin M → Fin r → ℝ := fun i j => (m.tbl i).θ j

/-- **Scope of the present model class.** `SpikedModelR.hθanti` makes every table order its
spikes strictly, and `hθnn` makes them nonnegative, so on the exactly aligned family the
paper's index is `ℓ_j = j` once component `j` carries a signal. The unordered case of
`alg:rank_r_stacksvd` (the worked example `θ_1 = (2,1)`, `θ_2 = (1,10)`) is therefore **not**
an `UnalignedModelR` with `R_i = 1`: it needs a permutation `R_i`, or a model without
`hθanti`. `Scalars.ellR` is kept general so that the statement does not have to change when
the model class is widened. The hypothesis `hex` is what the paper's `ℓ_j` needs: at least one
table carries component `j` (F8, 2026-09-05). -/
theorem ellR_thetaAligned (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) {j : Fin r} (hex : ∃ i, 0 < (m.tbl i).θ j) :
    Scalars.ellR m.thetaAligned c j = (j : ℕ) :=
  Scalars.ellR_eq_val hc (fun i k => (m.tbl i).hθnn k) (fun i => (m.tbl i).hθanti) hex

/-- `X_stack^(j)` (`eq:stacksvd_appXstack`): the stack weighted by the optimal weights of
component `j`. -/
noncomputable def stackXJ (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ)
    (j : Fin r) (N : ℕ) (ω : Ω N) : Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
  m.stackXW (Scalars.wStackR m.thetaAligned c j) N ω

/-- Gram matrix of `X_stack^(j)`. -/
noncomputable def stackGramJ (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ)
    (j : Fin r) (N : ℕ) (ω : Ω N) : Matrix (Fin (d N)) (Fin (d N)) ℝ :=
  (m.stackXJ c j N ω)ᵀ * m.stackXJ c j N ω

theorem isHermitian_stackGramJ (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (j : Fin r) (N : ℕ) (ω : Ω N) : (m.stackGramJ c j N ω).IsHermitian :=
  isHermitian_transpose_mul_self _

/-- `(v_jᵀ v̂_{j,stacksvd})²` in projector form: the squared overlap of `v_j` with the
`ℓ_j`-th right singular subspace of `X_stack^(j)`. Under simplicity of that eigenvalue it is
the paper's inner product squared (`stackOverlapJ_eq_inner_sq`). -/
noncomputable def stackOverlapJ (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (j : Fin r) (N : ℕ) (ω : Ω N) : ℝ :=
  overlapIdx (m.stackXJ c j N ω) (Scalars.ellR m.thetaAligned c j) (m.colVecG N j)

theorem stackOverlapJ_nonneg (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ)
    (j : Fin r) (N : ℕ) (ω : Ω N) : 0 ≤ m.stackOverlapJ c j N ω :=
  sq_nonneg _

/-- An orthogonal projector is a contraction and `‖v_j‖ = 1`, so the overlap is at most `1`.
Mirror: `overlap_le_norm_sq` (`Spectral.lean`). -/
theorem stackOverlapJ_le_one (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ)
    (j : Fin r) (N : ℕ) (ω : Ω N) : m.stackOverlapJ c j N ω ≤ 1 := by
  have h := overlapIdx_le_norm_sq (m.stackXJ c j N ω) (Scalars.ellR m.thetaAligned c j)
    (m.colVecG N j)
  rwa [m.norm_colVecG N j, one_pow] at h

/-- `v̂_{j,stacksvd}`: a unit eigenvector of `X_stack^(j)ᵀ X_stack^(j)` at the sorted index
`ℓ_j`, with the sign fixed by `⟪v̂_j, v_j⟫ ≥ 0`, as `vhatG` fixes it for one table. -/
noncomputable def vhatStackR (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ)
    (j : Fin r) (N : ℕ) (ω : Ω N) : EuclideanSpace ℝ (Fin (d N)) :=
  if 0 ≤ ⟪vEig (m.stackGramJ c j N ω) (m.isHermitian_stackGramJ c j N ω)
      (Scalars.ellR m.thetaAligned c j), m.colVecG N j⟫_ℝ then
    vEig (m.stackGramJ c j N ω) (m.isHermitian_stackGramJ c j N ω)
      (Scalars.ellR m.thetaAligned c j)
  else
    -vEig (m.stackGramJ c j N ω) (m.isHermitian_stackGramJ c j N ω)
      (Scalars.ellR m.thetaAligned c j)

/-- The sign flip does not change a squared inner product, so `vhatStackR` and the raw
eigenvector give the same overlaps. -/
theorem inner_vhatStackR_sq (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ)
    (j k : Fin r) (N : ℕ) (ω : Ω N) :
    ⟪m.vhatStackR c j N ω, m.colVecG N k⟫_ℝ ^ 2
      = ⟪vEig (m.stackGramJ c j N ω) (m.isHermitian_stackGramJ c j N ω)
          (Scalars.ellR m.thetaAligned c j), m.colVecG N k⟫_ℝ ^ 2 := by
  rw [vhatStackR]
  split_ifs with h
  · rfl
  · rw [inner_neg_left, neg_sq]

/-- The paper's sign convention `⟪v̂_j, v_j⟫ ≥ 0`. -/
theorem inner_vhatStackR_nonneg (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (j : Fin r) (N : ℕ) (ω : Ω N) :
    0 ≤ ⟪m.vhatStackR c j N ω, m.colVecG N j⟫_ℝ := by
  rw [vhatStackR]
  split_ifs with h
  · exact h
  · rw [inner_neg_left]
    linarith [not_le.mp h]

/-- `v̂_j` is a unit vector whenever the index `ℓ_j` is inside the dimension. -/
theorem norm_vhatStackR (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ)
    (j : Fin r) (N : ℕ) (ω : Ω N) (hk : Scalars.ellR m.thetaAligned c j < d N) :
    ‖m.vhatStackR c j N ω‖ = 1 := by
  have hkc : Scalars.ellR m.thetaAligned c j < Fintype.card (Fin (d N)) := by simpa using hk
  rw [vhatStackR]
  split_ifs with h
  · exact norm_vEig _ _ hkc
  · rw [norm_neg]
    exact norm_vEig _ _ hkc

/-- The projector form and the paper's inner-product form agree once the `ℓ_j`-th eigenvalue
of the Gram matrix is simple **at that one index**. Mirror: `thm_stacksvd_weighted_inner` at
rank 1. The hypothesis is `SimpleIdx` (`LinAlg/SpecIdxPerturb.lean`), not `SimpleSpec` at
`ℓ_j + 1`: the proof reads one index only, and the Gaussian discharge (Track E) supplies the
one-index form. -/
theorem stackOverlapJ_eq_inner_sq (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (j : Fin r) (N : ℕ) (ω : Ω N)
    (hs : SimpleIdx (m.stackGramJ c j N ω) (m.isHermitian_stackGramJ c j N ω)
      (Scalars.ellR m.thetaAligned c j)) :
    m.stackOverlapJ c j N ω = ⟪m.vhatStackR c j N ω, m.colVecG N j⟫_ℝ ^ 2 := by
  rw [inner_vhatStackR_sq, stackOverlapJ, overlapIdx]
  exact normSq_specProjIdx_eq_inner_sq hs _

/-- `‖Vᵀ V̂_stacksvd‖_F² = ∑_j ∑_k ⟪v̂_j, v_k⟫²`, the aggregate of the corollary
(`main_paper.tex:2344`), written on the entries so that no `d × r` matrix is built. -/
noncomputable def frobSqStackR (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (N : ℕ) (ω : Ω N) : ℝ :=
  ∑ j : Fin r, ∑ k : Fin r, ⟪m.vhatStackR c j N ω, m.colVecG N k⟫_ℝ ^ 2

theorem frobSqStackR_nonneg (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ)
    (N : ℕ) (ω : Ω N) : 0 ≤ m.frobSqStackR c N ω :=
  Finset.sum_nonneg fun _ _ => Finset.sum_nonneg fun _ _ => sq_nonneg _

/-- The diagonal of the Frobenius aggregate is the `r` per-component overlaps, once every
`ℓ_j`-th eigenvalue is simple at its own index (`SimpleIdx`). The off-diagonal terms are what
the limit theorem (item D2) must send to `0`. -/
theorem frobSqStackR_diag (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ)
    (N : ℕ) (ω : Ω N)
    (hs : ∀ j : Fin r, SimpleIdx (m.stackGramJ c j N ω) (m.isHermitian_stackGramJ c j N ω)
      (Scalars.ellR m.thetaAligned c j)) :
    m.frobSqStackR c N ω
      = (∑ j : Fin r, m.stackOverlapJ c j N ω)
        + ∑ j : Fin r, ∑ k ∈ Finset.univ.erase j,
            ⟪m.vhatStackR c j N ω, m.colVecG N k⟫_ℝ ^ 2 := by
  rw [frobSqStackR, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [m.stackOverlapJ_eq_inner_sq c j N ω (hs j), add_comm,
    Finset.sum_erase_add _ _ (Finset.mem_univ j)]

end AlignedStack

/-! ### 2. The hypothesis structure -/

/-- The rank-`r` twin of `MultiTableModel.HeteroLaw` (`StackSVDWeighted.lean:1267`), one
field per thing `thm:rank_r_stacksvd` reads. Every field is a statement about the `r`
weighted stacks `X_stack^(j)` of `eq:stacksvd_appXstack`.

`align` is the first display of the corollary in projector form. `crossProj` is the paper's
sentence "the columns of `V̂` will be asymptotically orthogonal" (`main_paper.tex:2318`),
which the second display needs and which `align` does not give. `simpleIdxJ` is what turns
the projector overlap into the paper's `⟪v̂_j, v_j⟫²`.

Two shape choices (D1 choice 11, Track E plan section 1.2, audit changes 3 and 4).

* `crossProj` is the **projector** form. It is stronger than the inner form by Bessel
  (`overlapIdx_ge_inner_sq`, `LinAlg/SpecIdx.lean:270`), so the inner form is a theorem
  (`HeteroLawR.cross`) that needs no simplicity. The reason for the shape is that the
  Gaussian discharge (Track E) proves projector forms; converting each instance through
  simplicity inside the discharge would be wasted work.
* `simpleIdxJ` asks simplicity at the **one** index `ℓ_j` (`SimpleIdx`,
  `LinAlg/SpecIdxPerturb.lean:484`), not `SimpleSpec` at `ℓ_j + 1`. The old field also ruled
  out a tie between two components above `θ̃_jj`, which no theorem here reads.

The paper's separation hypothesis `θ̃_jj ≠ θ̃_jk` for `k ≠ j` (`eq:stacksvd_apptildTheta`),
`0 < c_i` and `R_i = 1` are what make these three fields true; they are not read by any
theorem below, so they sit on a future Gaussian facade, not here. -/
structure HeteroLawR (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ) :
    Prop where
  /-- Per component: the squared overlap of `v_j` with the `ℓ_j`-th right singular subspace
  of `X_stack^(j)` tends to `γ_j` (`eq:stacksvd_gammak`, total by E3). -/
  align : ∀ j : Fin r, TendstoInProb μ (fun N ω => m.stackOverlapJ c j N ω)
    (Scalars.gammaR m.thetaAligned c j)
  /-- Asymptotic orthogonality, in projector form: the squared overlap of `v_k`, `k ≠ j`,
  with the `ℓ_j`-th right singular subspace of `X_stack^(j)` tends to `0`. -/
  crossProj : ∀ j k : Fin r, k ≠ j →
    TendstoInProb μ (fun N ω => overlapIdx (m.stackXJ c j N ω)
      (Scalars.ellR m.thetaAligned c j) (m.colVecG N k)) 0
  /-- The `ℓ_j`-th eigenvalue of the Gram matrix of `X_stack^(j)` is almost surely simple at
  that index. -/
  simpleIdxJ : ∀ (j : Fin r) (N : ℕ), ∀ᵐ ω ∂(μ N),
    SimpleIdx (m.stackGramJ c j N ω) (m.isHermitian_stackGramJ c j N ω)
      (Scalars.ellR m.thetaAligned c j)

end UnalignedModelR

end StackedSVD
