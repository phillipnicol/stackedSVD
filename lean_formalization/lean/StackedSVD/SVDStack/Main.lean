/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.SVDStack.Gram
import StackedSVD.RMT.Full
import StackedSVD.Scalars

/-!
# `thm:svd_stack_general`: the theorems

The last module of `StackedSVD.SVDStack`. It assembles `StackedSVD.SVDStack.Defs`,
`StackedSVD.SVDStack.Deterministic` and `StackedSVD.SVDStack.Gram` into the three statements
of the paper.

## Content

1. `thm_svd_stack_general`, the Layer 1 form
   `svdstackPerf → (βᵀ v_max(A_β))² / λ_max(A_β)`.
2. `thm_svd_stack_general_inner`, the paper's form `|⟨v̂_svdstack, v⟩|² → ...`.
3. `thm_svd_stack_general_zero`, the case `β_1 = 0`.
4. `thm_svd_stack_general_gaussian`, the Layer 2 corollary: the same conclusion from the
   proportional regime and the joint Gaussian law alone. `StackedSVD.RMT.Full` supplies
   `SpikedModel.singleTableLaw_of_gaussian` and `MultiTableModel.gaussianNoise_of_joint`
   supplies each table's marginal law. `RMT/*` never imports `SVDStack/*`, so the import adds
   no cycle.

The private helpers `svdstackPerf_le_sum` (the deterministic bound behind the zero case) and
`svdstackPerf_eq_of_goodEvent` (step (v) on the good event) are used only here.
`betaSq_nonneg` and `abs_norm_sq_topProj_sub_le` are public: `SVDStack/Weighted.lean` needs
both (`notes/archive/audit_weighted_2026-08-30.md`, section 3).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal Matrix.Norms.L2Operator

namespace StackedSVD

section SpectralHelpers

variable {d : ℕ}

/-- `⟪P w, w⟫ = ‖P w‖²` for an orthogonal projector. -/
private theorem inner_starProjection_self (K : Submodule ℝ (EuclideanSpace ℝ (Fin d)))
    (w : EuclideanSpace ℝ (Fin d)) : ⟪K.starProjection w, w⟫_ℝ = ‖K.starProjection w‖ ^ 2 := by
  have h0 : ⟪K.starProjection w, w - K.starProjection w⟫_ℝ = 0 :=
    K.sub_starProjection_mem_orthogonal w _ (K.starProjection_apply_mem w)
  have : ⟪K.starProjection w, w⟫_ℝ - ⟪K.starProjection w, K.starProjection w⟫_ℝ = 0 := by
    rw [← inner_sub_right]
    exact h0
  rw [← real_inner_self_eq_norm_sq]
  linarith

/-- `⟪P_top w, w⟫ = ‖P_top w‖²`, the `topProj` form of `inner_starProjection_self`. -/
private theorem inner_topProj_self (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (w : EuclideanSpace ℝ (Fin d)) : ⟪topProj A hA w, w⟫_ℝ = ‖topProj A hA w‖ ^ 2 :=
  inner_starProjection_self _ w

/-- `betaSq` is nonnegative for every `θ` and `c`. -/
theorem betaSq_nonneg (θ c : ℝ) : 0 ≤ betaSq θ c := by
  unfold betaSq
  split
  · apply div_nonneg
    · linarith
    · positivity
  · exact le_rfl

/-- The two squared norms of a projector applied to two vectors are close when the vectors
are: `|‖P y‖² - ‖P b‖²| ≤ ‖y - b‖ (‖y‖ + ‖b‖)`, by `‖P‖ ≤ 1`. -/
theorem abs_norm_sq_topProj_sub_le (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (y b : EuclideanSpace ℝ (Fin d)) :
    |‖topProj A hA y‖ ^ 2 - ‖topProj A hA b‖ ^ 2| ≤ ‖y - b‖ * (‖y‖ + ‖b‖) := by
  have hy : ‖topProj A hA y‖ ≤ ‖y‖ := Submodule.norm_starProjection_apply_le _ y
  have hb : ‖topProj A hA b‖ ≤ ‖b‖ := Submodule.norm_starProjection_apply_le _ b
  have hsub : ‖topProj A hA y - topProj A hA b‖ ≤ ‖y - b‖ := by
    rw [← map_sub]
    exact Submodule.norm_starProjection_apply_le _ _
  have h1 : |‖topProj A hA y‖ - ‖topProj A hA b‖| ≤ ‖y - b‖ :=
    (abs_norm_sub_norm_le _ _).trans hsub
  have hfac : ‖topProj A hA y‖ ^ 2 - ‖topProj A hA b‖ ^ 2
      = (‖topProj A hA y‖ - ‖topProj A hA b‖) * (‖topProj A hA y‖ + ‖topProj A hA b‖) := by
    ring
  have hnn : 0 ≤ ‖topProj A hA y‖ + ‖topProj A hA b‖ := by positivity
  rw [hfac, abs_mul, abs_of_nonneg hnn]
  exact mul_le_mul h1 (by linarith) hnn (by positivity)

end SpectralHelpers

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

section SVDStack

/-! ### The svdstack performance and `thm:svd_stack_general` -/

-- One table at least. `Fin M` must be nonempty: the shared `v` is read off table `0`.
variable [NeZero M]

/-- Deterministic bound behind the zero case. The top eigenvalue of `∑_i P_i` is at least `1`,
because `∑_i P_i ⪰ P_0` and `P_0` is a nonzero projector. Writing `u = P_top(∑_i P_i) v`, the
eigenvalue relation gives `λ ‖u‖² = ⟪u, ∑_i P_i v⟫ ≤ ‖u‖ ∑_i ‖P_i v‖`, so
`‖u‖² ≤ (∑_i ‖P_i v‖)²`. -/
private theorem svdstackPerf_le_sum (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) :
    m.svdstackPerf N ω ≤
      (∑ i, ‖topProj (m.tableGram i N ω) (m.isHermitian_tableGram i N ω) ((m.tbl 0).v N)‖) ^ 2 := by
  set v := (m.tbl 0).v N with hvdef
  set G := m.svdstackGram N ω with hGdef
  set hGh := m.isHermitian_svdstackGram N ω with hGhdef
  set Q : Fin M → EuclideanSpace ℝ (Fin (d N)) →L[ℝ] EuclideanSpace ℝ (Fin (d N)) :=
    fun i => topProj (m.tableGram i N ω) (m.isHermitian_tableGram i N ω) with hQdef
  set S : ℝ := ∑ i, ‖Q i v‖ with hSdef
  have hSnn : 0 ≤ S := Finset.sum_nonneg fun i _ => norm_nonneg _
  -- the operator of `∑_i P_i` is `∑_i Q i`
  have hop : ∀ x, toOp G x = ∑ i, Q i x := by
    have hmap : Matrix.toEuclideanLin G = ∑ i, ((Q i : _ →L[ℝ] _) : _ →ₗ[ℝ] _) := by
      rw [hGdef]
      change Matrix.toEuclideanLin (∑ i, m.P i N ω) = _
      rw [map_sum]
      refine Finset.sum_congr rfl fun i _ => ?_
      simp only [MultiTableModel.P, hQdef, LinearEquiv.apply_symm_apply]
    intro x
    rw [show toOp G x = Matrix.toEuclideanLin G x from rfl, hmap]
    simp
  set u := topProj G hGh v with hudef
  have humem : u ∈ topSpace G hGh := Submodule.starProjection_apply_mem _ v
  have hGu : toOp G u = lamMax G hGh • u := toOp_of_mem_topSpace humem
  have hsym : (toOp G).IsSymmetric := Matrix.isSymmetric_toEuclideanLin_iff.mpr hGh
  -- `λ ≥ 1`
  have hd0 : 0 < d N := (m.tbl 0).hd N
  have hone : (1 : ℝ) ≤ lamMax G hGh := by
    set w := vMax (m.tableGram 0 N ω) (m.isHermitian_tableGram 0 N ω) with hwdef
    have hwn : ‖w‖ = 1 := norm_vMax hd0 _ _
    have hQ0 : Q 0 w = w :=
      Submodule.starProjection_eq_self_iff.mpr (mem_topSpace_vMax hd0 _ _)
    have hsum : ⟪toOp G w, w⟫_ℝ = ∑ i, ‖Q i w‖ ^ 2 := by
      rw [hop w, sum_inner]
      exact Finset.sum_congr rfl fun i _ => inner_starProjection_self _ _
    have h1 : (1 : ℝ) ≤ ∑ i, ‖Q i w‖ ^ 2 := by
      have := Finset.single_le_sum (f := fun i => ‖Q i w‖ ^ 2)
        (fun i _ => sq_nonneg _) (Finset.mem_univ (0 : Fin M))
      rwa [hQ0, hwn, one_pow] at this
    have h2 := inner_toOp_self_le G hGh w
    rw [hsum, hwn] at h2
    nlinarith
  -- `λ ‖u‖² ≤ ‖u‖ S`
  have key : lamMax G hGh * ‖u‖ ^ 2 ≤ ‖u‖ * S := by
    have h1 : ⟪toOp G u, v⟫_ℝ = lamMax G hGh * ‖u‖ ^ 2 := by
      rw [hGu, real_inner_smul_left]
      congr 1
      exact inner_topProj_self G hGh v
    have h2 : ⟪toOp G u, v⟫_ℝ = ∑ i, ⟪u, Q i v⟫_ℝ := by
      rw [hsym u v, hop v, inner_sum]
    have h3 : ∑ i, ⟪u, Q i v⟫_ℝ ≤ ∑ i, ‖u‖ * ‖Q i v‖ :=
      Finset.sum_le_sum fun i _ => real_inner_le_norm _ _
    rw [← h1, h2]
    calc ∑ i, ⟪u, Q i v⟫_ℝ ≤ ∑ i, ‖u‖ * ‖Q i v‖ := h3
      _ = ‖u‖ * S := by rw [hSdef, Finset.mul_sum]
  -- conclude
  change ‖u‖ ^ 2 ≤ S ^ 2
  rcases eq_or_lt_of_le (norm_nonneg u) with hu0 | hu0
  · rw [← hu0]
    simpa using sq_nonneg S
  · have h4 : ‖u‖ * ‖u‖ ≤ ‖u‖ * S := by nlinarith
    have h5 : ‖u‖ ≤ S := le_of_mul_le_mul_left h4 hu0
    nlinarith

/-- Step (v) of the paper's proof on the good event: `svdstackPerf = ‖P_top(G_N) (Ṽ v)‖² /
λ_max(G_N)` with `G_N = Ṽ Ṽᵀ`. Route: `svdstackGram = Ṽᵀ Ṽ`, whose top eigenvalue is simple
and positive by transfer from `G_N` (`lamMax_mul_transpose_self_eq`,
`topSimple_transpose_mul_iff`); `svdstackPerf_eq_closed` gives the paper's ratio
`(xᵀ Ṽ v)² / (xᵀ G_N x)` with `x = v_max(G_N)`; the denominator is `λ_max(G_N)` and the
numerator is `‖P_top(G_N) (Ṽ v)‖²` by `norm_topProj_sq_eq_inner_sq`. -/
private theorem svdstackPerf_eq_of_goodEvent (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N)
    (h : m.goodEvent N ω) :
    m.svdstackPerf N ω =
      ‖topProj (m.Vt N ω * (m.Vt N ω)ᵀ) (isHermitian_mul_transpose_self (m.Vt N ω))
        (WithLp.toLp 2 ((m.Vt N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N))))‖ ^ 2 /
      lamMax (m.Vt N ω * (m.Vt N ω)ᵀ) (isHermitian_mul_transpose_self (m.Vt N ω)) := by
  obtain ⟨heq, hsimple, hpos⟩ := h
  set Vt := m.Vt N ω with hVt
  set v := (m.tbl 0).v N with hv
  have hA := isHermitian_transpose_mul_self Vt
  have hB := isHermitian_mul_transpose_self Vt
  have hlam : lamMax (Vtᵀ * Vt) hA = lamMax (Vt * Vtᵀ) hB :=
    lamMax_mul_transpose_self_eq Vt hpos
  have hposA : 0 < lamMax (Vtᵀ * Vt) hA := by
    rw [hlam]
    exact hpos
  have hsimpleA : TopSimple (Vtᵀ * Vt) hA :=
    (topSimple_transpose_mul_iff Vt hA hB hlam hpos).mpr hsimple
  have hperf : m.svdstackPerf N ω = ‖topProj (Vtᵀ * Vt) hA v‖ ^ 2 := by
    change ‖topProj (m.svdstackGram N ω) (m.isHermitian_svdstackGram N ω) v‖ ^ 2 = _
    rw [topProj_congr heq (m.isHermitian_svdstackGram N ω) hA]
  rw [hperf, svdstackPerf_eq_closed Vt v hsimpleA hposA, svdstackPerfClosed]
  have hM : 0 < M := NeZero.pos M
  set x := vMax (Vt * Vtᵀ) hB with hxdef
  have hxn : ‖x‖ = 1 := norm_vMax hM _ hB
  have hxmem : x ∈ topSpace (Vt * Vtᵀ) hB := mem_topSpace_vMax hM _ hB
  have hBx : toOp (Vt * Vtᵀ) x = lamMax (Vt * Vtᵀ) hB • x := toOp_of_mem_topSpace hxmem
  have hden : WithLp.ofLp x ⬝ᵥ (Vt * Vtᵀ).mulVec (WithLp.ofLp x) = lamMax (Vt * Vtᵀ) hB := by
    have : WithLp.ofLp x ⬝ᵥ (Vt * Vtᵀ).mulVec (WithLp.ofLp x)
        = ⟪x, toOp (Vt * Vtᵀ) x⟫_ℝ := by
      rw [real_inner_eq_dotProduct]
      rfl
    rw [this, hBx, real_inner_smul_right, real_inner_self_eq_norm_sq, hxn]
    ring
  have hnum : WithLp.ofLp x ⬝ᵥ Vt.mulVec (WithLp.ofLp v)
      = ⟪x, (WithLp.toLp 2 (Vt.mulVec (WithLp.ofLp v)) : EuclideanSpace ℝ (Fin M))⟫_ℝ := by
    rw [real_inner_eq_dotProduct]
  rw [hden, hnum, norm_topProj_sq_eq_inner_sq hsimple hxmem hxn]

omit [NeZero M] in
/-- `thm:svd_stack_general`, Layer 1 form (`law` and `hI` in the signature: the single-table
law structure for each table, and the independence of the tables; since 2026-09-02 (L4) the
noise enters only through `IndepNoise`, the Gaussian law only through `law`). `hthr` is
the paper's simplicity condition `β_2 > 0` in unsorted form: two distinct tables have a
positive `β`. It gives the spectral gap `λ_1(A_β) - λ_2(A_β) ≥ β_i β_j > 0` through
`abeta_gap`, and it forces `2 ≤ M`. No ordering of `β` is assumed, so the theorem covers an
unsorted model with no permutation lemma. `hc` gives `0 ≤ β_i < 1` through the formula for
`betaSq`, so the earlier hypothesis `hβ01` is not needed. -/
theorem thm_svd_stack_general [NeZero M] (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ i j, i ≠ j ∧ 0 < β i ∧ 0 < β j)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) (svdstackLimit β) := by
  obtain ⟨i, j, hij, hi, hj⟩ := id hthr
  have h1 : 1 < Fintype.card (Fin M) := one_lt_card_fin_of_ne hij
  have hβ01 : ∀ k, 0 ≤ β k ∧ β k < 1 := fun k => by
    rw [hβdef k]
    exact beta_mem_Ico (hc k)
  have hA : (Abeta β).IsHermitian := isHermitian_Abeta β
  have hgap := abeta_gap β hβ01 hij hi hj
  have hγ : 0 < lamMax (Abeta β) hA - hA.eigenvalues₀ ⟨1, h1⟩ :=
    lt_of_lt_of_le (mul_pos hi hj) hgap
  have hsimple : TopSimple (Abeta β) hA := topSimple_of_gap hA h1 hγ
  have hApos : 0 < lamMax (Abeta β) hA := lt_of_lt_of_le zero_lt_one (one_le_lamMax_Abeta β i)
  set bvec : EuclideanSpace ℝ (Fin M) := WithLp.toLp 2 β with hbvec
  have hsymm : ∀ (N : ℕ) (ω : Ω N), (m.Vt N ω * (m.Vt N ω)ᵀ).IsHermitian :=
    fun N ω => isHermitian_mul_transpose_self (m.Vt N ω)
  -- (a) `λ_max(G_N) → λ_max(A_β)`, by Weyl
  have hlam : TendstoInProb μ (fun N ω =>
      lamMax (m.Vt N ω * (m.Vt N ω)ᵀ) (isHermitian_mul_transpose_self (m.Vt N ω)))
      (lamMax (Abeta β) hA) :=
    lamMax_tendstoInProb hA hsymm (m.gramEntries c β hβdef law hI)
  -- (b) `‖P_top(G_N) β‖² → ⟪v_max(A_β), β⟫²`, by `lem_entrywise_conv_eigenvec`
  have hnumβ : TendstoInProb μ (fun N ω =>
      ‖topProj (m.Vt N ω * (m.Vt N ω)ᵀ) (isHermitian_mul_transpose_self (m.Vt N ω)) bvec‖ ^ 2)
      (⟪vMax (Abeta β) hA, bvec⟫_ℝ ^ 2) :=
    lem_entrywise_conv_eigenvec β (fun N ω => m.Vt N ω * (m.Vt N ω)ᵀ) hsymm
      (m.gramEntries c β hβdef law hI) hsimple bvec
  -- (c) `y_N := Ṽ v → β`, so `‖P_top(G_N) y_N‖² - ‖P_top(G_N) β‖² → 0`
  have hy : TendstoInProbPi μ (fun N ω => (m.Vt N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N))) β := by
    intro k
    rw [hβdef k]
    exact m.align c law k
  have hyβ : TendstoInProb μ (fun N ω =>
      ‖(WithLp.toLp 2 ((m.Vt N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N))) :
          EuclideanSpace ℝ (Fin M)) - bvec‖ *
        (‖(WithLp.toLp 2 ((m.Vt N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N))) :
          EuclideanSpace ℝ (Fin M))‖ + ‖bvec‖)) 0 := by
    have hφ : Continuous (fun w : Fin M → ℝ =>
        ‖(WithLp.toLp 2 w : EuclideanSpace ℝ (Fin M)) - bvec‖ *
          (‖(WithLp.toLp 2 w : EuclideanSpace ℝ (Fin M))‖ + ‖bvec‖)) := by
      fun_prop
    have h := TendstoInProbPi.comp_continuous hφ.continuousAt hy
    rw [hbvec, sub_self, norm_zero, zero_mul] at h
    exact h
  have hdiff : TendstoInProb μ (fun N ω =>
      ‖topProj (m.Vt N ω * (m.Vt N ω)ᵀ) (isHermitian_mul_transpose_self (m.Vt N ω))
        (WithLp.toLp 2 ((m.Vt N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N))))‖ ^ 2 -
      ‖topProj (m.Vt N ω * (m.Vt N ω)ᵀ) (isHermitian_mul_transpose_self (m.Vt N ω)) bvec‖ ^ 2)
      0 := by
    refine TendstoInProb.of_le (fun N => ?_) hyβ
    filter_upwards with ω
    rw [sub_zero]
    exact abs_norm_sq_topProj_sub_le _ _ _ _
  have hnum : TendstoInProb μ (fun N ω =>
      ‖topProj (m.Vt N ω * (m.Vt N ω)ᵀ) (isHermitian_mul_transpose_self (m.Vt N ω))
        (WithLp.toLp 2 ((m.Vt N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N))))‖ ^ 2)
      (⟪vMax (Abeta β) hA, bvec⟫_ℝ ^ 2) := by
    have h := hdiff.add hnumβ
    rw [zero_add] at h
    refine h.congr fun N => ?_
    filter_upwards with ω
    ring
  -- (d) the quotient, and the identification of the limit
  have hg := hnum.div hlam hApos.ne'
  have hlimit : svdstackLimit β = ⟪vMax (Abeta β) hA, bvec⟫_ℝ ^ 2 / lamMax (Abeta β) hA := by
    rw [svdstackLimit, real_inner_eq_dotProduct, dotProduct_comm]
  rw [hlimit]
  -- (e) transfer from the good event, where `svdstackPerf` equals the quotient
  refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_ hg
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds
    (m.tendsto_measure_not_goodEvent c β hc hβdef hthr law hI) (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  change ¬ m.goodEvent N ω
  intro hgood
  exact hω (m.svdstackPerf_eq_of_goodEvent N ω hgood)

/-- `thm:svd_stack_general` in the paper's form `|⟨v̂_svdstack, v⟩|² → ...`, for any selection
`vhat` of a unit top eigenvector of `∑_i P_i` (for instance `svdstackEst`, or its negative).
No measurability of `vhat` is needed: `TendstoInProb` is a statement about the outer measure
of arbitrary sets, and the transfer from `svdstackPerf` goes through the event
`{¬ TopSimple (svdstackGram N ω)}` of `topSimple_svdstackGram`
(`TendstoInProb.of_tendsto_measure_ne_of_tendsto`), outside of which the two agree by
`norm_topProj_sq_eq_inner_sq`. -/
theorem thm_svd_stack_general_inner (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ i j, i ≠ j ∧ 0 < β i ∧ 0 < β j)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise)
    (vhat : (N : ℕ) → Ω N → EuclideanSpace ℝ (Fin (d N)))
    (hmem : ∀ N ω, vhat N ω ∈ topSpace (m.svdstackGram N ω) (m.isHermitian_svdstackGram N ω))
    (hnorm : ∀ N ω, ‖vhat N ω‖ = 1) :
    TendstoInProb μ (fun N ω => ⟪vhat N ω, (m.tbl 0).v N⟫_ℝ ^ 2) (svdstackLimit β) := by
  refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_
    (m.thm_svd_stack_general c β hc hβdef hthr law hI)
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds
    (m.topSimple_svdstackGram c β hc hβdef hthr law hI) (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  have hω' : ⟪vhat N ω, (m.tbl 0).v N⟫_ℝ ^ 2 ≠ m.svdstackPerf N ω := hω
  change ¬ TopSimple (m.svdstackGram N ω) (m.isHermitian_svdstackGram N ω)
  intro hsimple
  exact hω' (norm_topProj_sq_eq_inner_sq hsimple (hmem N ω) (hnorm N ω) _).symm

-- `hI` is not used by the proof below: the bound `svdstackPerf ≤ (∑_i ‖P_i v‖)²` is
-- deterministic, and every term of the sum tends to `0` by `SingleTableLaw.align` alone. The
-- hypothesis stays in the statement to match `thm_svd_stack_general`.
set_option linter.unusedVariables false in
/-- `thm:svd_stack_general`, the case `β_1 = 0`. Under the paper's ordering that means every
`β_i = 0`, which is how the unsorted statement writes it. The paper leaves `β_1 > 0 = β_2`
open (`sec:svdstack_threshold`). -/
theorem thm_svd_stack_general_zero (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)) (hβ0 : ∀ i, β i = 0)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) 0 := by
  have hb : ∀ i, betaSq (m.tbl i).θ (c i) = 0 := by
    intro i
    have h : (0 : ℝ) = Real.sqrt (betaSq (m.tbl i).θ (c i)) := by
      rw [← hβ0 i, hβdef i, beta]
    have hle : betaSq (m.tbl i).θ (c i) ≤ 0 := Real.sqrt_eq_zero'.mp h.symm
    linarith [betaSq_nonneg (m.tbl i).θ (c i)]
  set F : (N : ℕ) → Ω N → Fin M → ℝ := fun N ω i =>
    Real.sqrt (overlap ((m.tbl i).X N ω) ((m.tbl 0).v N)) with hFdef
  have hFconv : TendstoInProbPi μ F fun _ => (0 : ℝ) := by
    intro i
    have h2 : TendstoInProb μ (fun N ω => overlap ((m.tbl i).X N ω) ((m.tbl 0).v N))
        (betaSq (m.tbl i).θ (c i)) := by
      have := (law i).align
      simpa [m.hv i 0] using this
    rw [hb i] at h2
    have h3 := h2.comp_continuous (φ := Real.sqrt) Real.continuous_sqrt.continuousAt
    simpa [hFdef] using h3
  have hlim : TendstoInProb μ (fun N ω => (∑ i, F N ω i) ^ 2) 0 := by
    have := TendstoInProbPi.comp_continuous (φ := fun w : Fin M → ℝ => (∑ i, w i) ^ 2)
      (Continuous.continuousAt (by fun_prop)) hFconv
    simpa using this
  refine TendstoInProb.of_le (g := fun N ω => (∑ i, F N ω i) ^ 2) (fun N => ?_) hlim
  filter_upwards with ω
  have hle := m.svdstackPerf_le_sum N ω
  have hEq : ∀ i : Fin M, F N ω i =
      ‖topProj (m.tableGram i N ω) (m.isHermitian_tableGram i N ω) ((m.tbl 0).v N)‖ := by
    intro i
    rw [hFdef]
    exact Real.sqrt_sq (norm_nonneg _)
  have hnn : (0 : ℝ) ≤ m.svdstackPerf N ω := by
    change (0 : ℝ) ≤
      ‖topProj (m.svdstackGram N ω) (m.isHermitian_svdstackGram N ω) ((m.tbl 0).v N)‖ ^ 2
    positivity
  rw [sub_zero, abs_of_nonneg hnn]
  calc m.svdstackPerf N ω
      ≤ (∑ i, ‖topProj (m.tableGram i N ω) (m.isHermitian_tableGram i N ω)
          ((m.tbl 0).v N)‖) ^ 2 := hle
    _ = (∑ i, F N ω i) ^ 2 := by
        congr 1
        exact Finset.sum_congr rfl fun i _ => (hEq i).symm

omit [NeZero M] in
/-- `thm:svd_stack_general` (`main_paper.tex:479`) with every random matrix theory hypothesis
discharged. The tables are Gaussian and independent (`hG`) and each is in the proportional
regime (`hreg`); no `SingleTableLaw` is assumed. The proof applies
`SpikedModel.singleTableLaw_of_gaussian` to table `i`, whose Gaussian marginal comes from
`gaussianNoise_of_joint`. This is the paper's unweighted svdstack theorem for Gaussian noise,
with nothing left as a hypothesis except the model itself. -/
theorem thm_svd_stack_general_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ i j, i ≠ j ∧ 0 < β i ∧ 0 < β j)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) (svdstackLimit β) :=
  m.thm_svd_stack_general c β hc hβdef hthr
    (fun i => SpikedModel.singleTableLaw_of_gaussian (hc i) (m.tbl i) (hreg i)
      (m.gaussianNoise_of_joint hG i)) hG.indepNoise

omit [NeZero M] in
/-- `thm:svd_stack_general`, the case `β_1 = 0` (`sec:svdstack_threshold`), with every random
matrix theory hypothesis discharged: the tables are Gaussian and independent (`hG`) and each is
in the proportional regime (`hreg`); no `SingleTableLaw` is assumed. The wrapper adds
`hc : ∀ i, 0 < c i`, which the base theorem does not need but which
`SpikedModel.singleTableLaw_of_gaussian` does. -/
theorem thm_svd_stack_general_zero_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)) (hβ0 : ∀ i, β i = 0)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) 0 :=
  m.thm_svd_stack_general_zero c β hβdef hβ0
    (fun i => SpikedModel.singleTableLaw_of_gaussian (hc i) (m.tbl i) (hreg i)
      (m.gaussianNoise_of_joint hG i)) hG.indepNoise

/-- `thm:svd_stack_general` in the paper's form `|⟨v̂_svdstack, v⟩|² → ...`, with every random
matrix theory hypothesis discharged: the tables are Gaussian and independent (`hG`) and each is
in the proportional regime (`hreg`); no `SingleTableLaw` is assumed. -/
theorem thm_svd_stack_general_inner_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ i j, i ≠ j ∧ 0 < β i ∧ 0 < β j)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (vhat : (N : ℕ) → Ω N → EuclideanSpace ℝ (Fin (d N)))
    (hmem : ∀ N ω, vhat N ω ∈ topSpace (m.svdstackGram N ω) (m.isHermitian_svdstackGram N ω))
    (hnorm : ∀ N ω, ‖vhat N ω‖ = 1) :
    TendstoInProb μ (fun N ω => ⟪vhat N ω, (m.tbl 0).v N⟫_ℝ ^ 2) (svdstackLimit β) :=
  m.thm_svd_stack_general_inner c β hc hβdef hthr
    (fun i => SpikedModel.singleTableLaw_of_gaussian (hc i) (m.tbl i) (hreg i)
      (m.gaussianNoise_of_joint hG i)) hG.indepNoise vhat hmem hnorm

/-- `thm:simple_thm1`, svdstack half, with every random matrix theory hypothesis discharged.
The tables share signal `θ₀` and regime `c₀`, above the detectability threshold `c₀ < θ₀⁴`
(`hdet`), with at least two tables (`hM2`); the limit is the paper's closed form
`1 - (c₀+θ₀²)/(Mθ₀⁴+θ₀²-(M-1)c₀)` (`Scalars.simple_thm1_svdstack`). The two distinct indices
needed for `hthr` come from `hM2`, and `beta_pos_of_thr` gives each table's `β` positive. -/
theorem thm_simple_thm1_svdstack_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (θ₀ c₀ : ℝ) (hc : 0 < c₀) (hdet : c₀ < θ₀ ^ 4)
    (hM2 : 2 ≤ M) (hθ : ∀ i, (m.tbl i).θ = θ₀)
    (hreg : ∀ i, (m.tbl i).Regime c₀) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω)
      (1 - (c₀ + θ₀ ^ 2) / ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2 - ((M : ℝ) - 1) * c₀)) := by
  have hM : 0 < M := by omega
  have hij : (⟨0, by omega⟩ : Fin M) ≠ (⟨1, by omega⟩ : Fin M) := by
    intro h
    have := congrArg Fin.val h
    simp at this
  -- `beta_pos_of_thr` lives in `SVDStack/Weighted.lean`, which imports this file, so its
  -- short proof is inlined here instead (deviation, see `notes/archive/agent_reports/facade_S.md`).
  have hβpos : 0 < beta θ₀ c₀ := by
    have hθ2 : 0 < θ₀ ^ 2 := by nlinarith [sq_nonneg θ₀, sq_nonneg (θ₀ ^ 2)]
    have hden : (0 : ℝ) < θ₀ ^ 4 + θ₀ ^ 2 := by nlinarith
    rw [beta, Real.sqrt_pos, betaSq, if_pos hdet]
    exact div_pos (by linarith) hden
  have hthr : ∃ i j : Fin M, i ≠ j ∧ 0 < beta θ₀ c₀ ∧ 0 < beta θ₀ c₀ :=
    ⟨⟨0, by omega⟩, ⟨1, by omega⟩, hij, hβpos, hβpos⟩
  have h := m.thm_svd_stack_general_gaussian (fun _ => c₀) (fun _ => beta θ₀ c₀)
    (fun _ => hc) (fun i => by rw [hθ i]) hthr hreg hG
  rwa [Scalars.simple_thm1_svdstack hM hc, if_pos hdet] at h

end SVDStack

end MultiTableModel

end StackedSVD
