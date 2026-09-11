/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.SVDStack.Deterministic
import StackedSVD.SVDStack.DelocDir

/-!
# `thm:svd_stack_general`: the consequences of the single-table laws

The third module of `StackedSVD.SVDStack`. Every result here reads the model
`MultiTableModel` and the hypotheses `SingleTableLaw` and `IndepNoise` (independent tables
with arbitrary laws; `JointGaussianNoise` implies it, `JointGaussianNoise.indepNoise`).

## Content

1. `svdstackGram_eq`: `∑_i P_i = Ṽᵀ Ṽ`, almost surely for each `N`.
2. The `delocDir` machinery of `lem:delocalization`: `perpOf`, `perpTopVec`, `delocDir`,
   `pi_pair_le` and `measure_deloc_le`. `delocDir v X` is a measurable unit vector orthogonal
   to `v` that is parallel to `(I - v vᵀ) v̂(X)` whenever the top eigenvalue of `Xᵀ X` is
   simple. It makes the random test direction of `lem:delocalization` a measurable function
   of one table alone, so that Fubini on the product law reduces the cross-table term to
   `delocUniform` at a fixed direction.
3. `align`, `lem_delocalization` and `gram`, the signed limits of `⟪v̂_i, v⟫` and
   `⟪v̂_i, v̂_j⟫`.
3a. `hasLaw_eval_of_hasLaw_pi` and `gaussianNoise_of_joint`: the marginal of the product law
   at one table, which every Layer 2 corollary of this file needs.
4. The good event `goodEvent` (`∑_i P_i = Ṽᵀ Ṽ`, and `Ṽ Ṽᵀ` has a simple positive top
   eigenvalue), `tendsto_measure_not_goodEvent`, and `topSimple_svdstackGram`.

The `delocDir` machinery, `measure_deloc_le` and `align_inner` are public, not `private`:
the unaligned rank-`r` results of `StackedSVD.RankR` run the same argument on `UnalignedModel`
and Lean 4 `private` is module scoped. Only `topSimple_svdstackGram_of_goodEvent` stays
`private`. The congruence helper is `topSimple_congr_iff` in `Spectral.lean`.

All limit statements about `Ṽ` are signed, not squared. A squared `gram` does not determine
the conclusion: two limit matrices with the same squared entries give svdstack limits
`0.529412` and `0` (`notes/archive/audit_scope_2026-08-29.md` section 6.1). The sign convention of
`vhat` is what makes the signed form correct.

The `delocDir`, `perpOf`, `perpTopVec`, `pi_pair_le` and `measure_deloc_le` machinery of
item 2 above, and the `PiMarginal` and `TableFamily` sections, moved to
`SVDStack/DelocDir.lean` (F28, 2026-09-08).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal Matrix.Norms.L2Operator


namespace StackedSVD


namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

section SVDStack

/-! ### `svdstackGram` and `Ṽᵀ Ṽ` -/

/-- `P_i = v̂_i v̂_iᵀ` whenever the top eigenvalue of table `i` is simple. Extracted from the
proof of `svdstackGram_eq` so that the weighted twin `svdstackGramW_eq`
(`SVDStack/Weighted.lean`) can reuse it (`notes/archive/audit_weighted_2026-08-30.md`, section 3).
-/
theorem P_eq_vecMulVec (m : MultiTableModel μ M n d) (i : Fin M) (N : ℕ) (ω : Ω N)
    (hsimple : TopSimple (m.tableGram i N ω) (m.isHermitian_tableGram i N ω)) :
    m.P i N ω = Matrix.vecMulVec (WithLp.ofLp (m.vhat i N ω)) (WithLp.ofLp (m.vhat i N ω)) := by
  have hd0 : 0 < d N := (m.tbl i).hd N
  have hen : ‖m.vhat i N ω‖ = 1 := by
    rw [MultiTableModel.vhat]
    split
    · exact norm_vMax hd0 _ _
    · rw [norm_neg]; exact norm_vMax hd0 _ _
  have hmem : m.vhat i N ω ∈ topSpace (m.tableGram i N ω) (m.isHermitian_tableGram i N ω) := by
    rw [MultiTableModel.vhat]
    split
    · exact mem_topSpace_vMax hd0 _ _
    · exact Submodule.neg_mem _ (mem_topSpace_vMax hd0 _ _)
  have hpr := topProj_eq_rankOne hsimple hmem hen
  have hrank : topProj (m.tableGram i N ω) (m.isHermitian_tableGram i N ω)
      = InnerProductSpace.rankOne ℝ (m.vhat i N ω) (m.vhat i N ω) := by
    ext w
    rw [hpr w]
    simp
  rw [MultiTableModel.P, hrank, InnerProductSpace.symm_toEuclideanLin_rankOne]
  simp

/-- `∑_i P_i = Ṽᵀ Ṽ` on the almost sure event that every table has a simple top eigenvalue.
There `P_i = v̂_i v̂_iᵀ`, and the sum of the rank-one terms is `Ṽᵀ Ṽ`. The event is a.s. for
each `N` by `SingleTableLaw.topSimple`, so the statement is a.e. per `N`, not a limit. -/
theorem svdstackGram_eq (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (N : ℕ) :
    ∀ᵐ ω ∂(μ N), m.svdstackGram N ω = (m.Vt N ω)ᵀ * m.Vt N ω := by
  have hall : ∀ᵐ ω ∂(μ N), ∀ i : Fin M,
      TopSimple (m.tableGram i N ω) (m.isHermitian_tableGram i N ω) :=
    ae_all_iff.mpr fun i => (law i).topSimple N
  filter_upwards [hall] with ω hω
  have hP : ∀ i : Fin M, m.P i N ω =
      Matrix.vecMulVec (WithLp.ofLp (m.vhat i N ω)) (WithLp.ofLp (m.vhat i N ω)) :=
    fun i => m.P_eq_vecMulVec i N ω (hω i)
  ext k l
  rw [MultiTableModel.svdstackGram, Matrix.sum_apply, Matrix.mul_apply]
  exact Finset.sum_congr rfl fun i _ => by rw [hP i, Matrix.vecMulVec_apply]; rfl

/-! ### Basic facts about `vhat` -/

/-- `v̂_i` is a unit vector. -/
theorem norm_vhat (m : MultiTableModel μ M n d) (i : Fin M) (N : ℕ) (ω : Ω N) :
    ‖m.vhat i N ω‖ = 1 := by
  rw [MultiTableModel.vhat]
  split_ifs with h
  · exact norm_vMax ((m.tbl i).hd N) _ _
  · rw [norm_neg]
    exact norm_vMax ((m.tbl i).hd N) _ _

/-- `v̂_i` lies in the top eigenspace of `X_iᵀ X_i`. -/
theorem mem_topSpace_vhat (m : MultiTableModel μ M n d) (i : Fin M) (N : ℕ) (ω : Ω N) :
    m.vhat i N ω ∈ topSpace (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω)
      (isHermitian_transpose_mul_self ((m.tbl i).X N ω)) := by
  rw [MultiTableModel.vhat]
  split_ifs with h
  · exact mem_topSpace_vMax ((m.tbl i).hd N) _ _
  · exact Submodule.neg_mem _ (mem_topSpace_vMax ((m.tbl i).hd N) _ _)

/-- The sign convention of `vhat`: `⟪v̂_i, v⟫ ≥ 0`. -/
theorem inner_vhat_nonneg (m : MultiTableModel μ M n d) (i : Fin M) (N : ℕ) (ω : Ω N) :
    0 ≤ ⟪m.vhat i N ω, (m.tbl i).v N⟫_ℝ := by
  rw [MultiTableModel.vhat]
  split_ifs with h
  · exact h
  · have hneg : ⟪-vMax (m.tableGram i N ω) (m.isHermitian_tableGram i N ω),
        (m.tbl i).v N⟫_ℝ
        = -⟪vMax (m.tableGram i N ω) (m.isHermitian_tableGram i N ω), (m.tbl i).v N⟫_ℝ := by
      simp
    rw [hneg]
    linarith [not_le.mp h]

/-! ### The noise-to-data map -/

/-- Table `k` as a function of its noise matrix: `X_k = θ_k u_k v^T + d^{-1/2} Z`. The table
form of `SpikedModel.dataOf`. -/
noncomputable def dataOf (m : MultiTableModel μ M n d) (k : Fin M) (N : ℕ)
    (Z : Matrix (Fin (n k N)) (Fin (d N)) ℝ) : Matrix (Fin (n k N)) (Fin (d N)) ℝ :=
  (m.tbl k).dataOf N Z

/-- `dataOf` at the model noise is the data matrix of that table. -/
theorem dataOf_eq (m : MultiTableModel μ M n d) (k : Fin M) (N : ℕ) (ω : Ω N) :
    m.dataOf k N ((m.tbl k).Z N ω) = (m.tbl k).X N ω := rfl

/-- The noise-to-data map of one table is measurable. -/
theorem measurable_dataOf (m : MultiTableModel μ M n d) (k : Fin M) (N : ℕ) :
    Measurable (m.dataOf k N) := (m.tbl k).measurable_dataOf N

/-- `gaussianMatrix p q` is a probability measure. The proof is the instance
`instIsProbabilityMeasureGaussianMatrix` of `Prob/GaussianMatrix.lean`; this file kept its own
copy until cleanup wave 1 (2026-08-30), which imported that file instead. The name stays
because `RankR/Defs.lean` and `ThetaEst.lean` call it; a later wave that may edit those two
deletes it and uses the instance. -/
theorem isProbabilityMeasure_gaussianMatrix (p q : ℕ) :
    IsProbabilityMeasure (gaussianMatrix p q) :=
  instIsProbabilityMeasureGaussianMatrix p q

/-- Marginal of the joint law: `JointGaussianNoise` gives every table its own
`SpikedModel.GaussianNoise`. This is what lets a Layer 2 corollary discharge
`SingleTableLaw` table by table with `SpikedModel.singleTableLaw_of_gaussian`. -/
theorem gaussianNoise_of_joint (m : MultiTableModel μ M n d) (hG : m.JointGaussianNoise)
    (i : Fin M) : (m.tbl i).GaussianNoise :=
  gaussianNoise_of_joint_pi m.tbl hG i

/-! ### `IndepNoise` and the Fubini step of `lem:delocalization` -/

/-- The paper's `assum:general_noise` as far as the rank 1 SVDstack theorems read it: the
tables are independent, each with some probability law on its noise matrix.
`JointGaussianNoise` is the Gaussian case. Every Layer 1 theorem of the SVDstack family
(`lem_delocalization`, `thm_svd_stack_general`, `thm_svdstack_weighted`, and the bound of
`SVDStack/Rayleigh.lean`) reads the noise only through this predicate, because the single
place that needs it is the Fubini step `measure_deloc_le`. `thm:theta_est` (`ThetaEst.lean`)
is not in this family: its item P uses the Gaussian second moments. Same shape as
`UnalignedModelR.IndepNoise` (`RankR/GramR.lean`); added 2026-09-02 (L4). -/
def IndepNoise (m : MultiTableModel μ M n d) : Prop :=
  ∃ ν : (N : ℕ) → (i : Fin M) → Measure (Matrix (Fin (n i N)) (Fin (d N)) ℝ),
    (∀ N i, IsProbabilityMeasure (ν N i)) ∧
      ∀ N, HasLaw (fun ω i => (m.tbl i).Z N ω) (Measure.pi (ν N)) (μ N)

/-- Independent Gaussian tables are independent tables. -/
theorem JointGaussianNoise.indepNoise {m : MultiTableModel μ M n d}
    (hG : m.JointGaussianNoise) : m.IndepNoise :=
  ⟨fun N i => gaussianMatrix (n i N) (d N), fun _ _ => inferInstance, hG⟩

/-- The Fubini step of `lem:delocalization`, on a `MultiTableModel`. The whole argument is
`measure_deloc_le_of_pi_indep`, which reads only the family of tables and the product law. -/
theorem measure_deloc_le (m : MultiTableModel μ M n d) (hI : m.IndepNoise)
    {i j : Fin M} (hij : i ≠ j) (N : ℕ) {η : ℝ} (hη : 0 < η) :
    μ N {ω | η ≤ overlap ((m.tbl i).X N ω)
        (delocDir ((m.tbl i).v N) ((m.tbl j).X N ω))}
      ≤ ⨆ w ∈ (m.tbl i).orthUnit N, μ N {ω | η ≤ overlap ((m.tbl i).X N ω) w} := by
  obtain ⟨ν, hν, hind⟩ := hI
  exact measure_deloc_le_of_pi_indep m.tbl (hν := hν) hind hij N hη


/-! ### `align`, `lem:delocalization` and `gram` -/

/-- `⟪v̂_i, v_i⟫ → β_i`, the signed form of `SingleTableLaw.align`. The sign convention of
`vhat` turns the square root of the overlap into the inner product itself. -/
theorem align_inner (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (i : Fin M) :
    TendstoInProb μ (fun N ω => ⟪m.vhat i N ω, (m.tbl i).v N⟫_ℝ)
      (beta (m.tbl i).θ (c i)) := by
  have h1 : TendstoInProb μ
      (fun N ω => Real.sqrt (overlap ((m.tbl i).X N ω) ((m.tbl i).v N)))
      (Real.sqrt (betaSq (m.tbl i).θ (c i))) :=
    ((law i).align).comp_continuous (φ := Real.sqrt) Real.continuous_sqrt.continuousAt
  refine h1.congr fun N => ?_
  filter_upwards [(law i).topSimple N] with ω hω
  rw [overlap_eq_inner_sq _ _ hω (m.mem_topSpace_vhat i N ω) (m.norm_vhat i N ω),
    Real.sqrt_sq_eq_abs, abs_of_nonneg (m.inner_vhat_nonneg i N ω)]

/-! ### Consequences of the single-table laws -/

/-- `lem:delocalization` (`main_paper.tex:1111`), signed: for `i ≠ j`,
`⟪v̂_i, v̂_j⟫ → β_i β_j`. The paper splits
`v̂_iᵀ v̂_j = v̂_iᵀ v vᵀ v̂_j + v̂_iᵀ (I - v vᵀ) v̂_j` (`main_paper.tex:1123`). The first term
tends to `β_i β_j` by `prop:single_table` part 1 and the sign convention of `vhat`; the second
tends to `0` by `delocUniform` of table `i` applied to the random direction of table `j`,
which needs Fubini over the product law, hence `hI`. -/
theorem lem_delocalization (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise)
    {i j : Fin M} (hij : i ≠ j) :
    TendstoInProb μ (fun N ω => ⟪m.vhat i N ω, m.vhat j N ω⟫_ℝ)
      (beta (m.tbl i).θ (c i) * beta (m.tbl j).θ (c j)) := by
  have hAj : TendstoInProb μ (fun N ω => ⟪m.vhat j N ω, (m.tbl i).v N⟫_ℝ)
      (beta (m.tbl j).θ (c j)) := by
    refine (m.align_inner c law j).congr fun N => ?_
    filter_upwards with ω
    rw [m.hv j i N]
  have hA := (m.align_inner c law i).mul hAj
  have hB : TendstoInProb μ
      (fun N ω => ⟪m.vhat i N ω, perpOf ((m.tbl i).v N) (m.vhat j N ω)⟫_ℝ) 0 := by
    refine TendstoInProb.of_le
      (g := fun N ω => Real.sqrt (overlap ((m.tbl i).X N ω)
        (delocDir ((m.tbl i).v N) ((m.tbl j).X N ω)))) (fun N => ?_) ?_
    · filter_upwards [(law j).topSimple N] with ω hω
      rw [sub_zero]
      have h1 : |⟪m.vhat i N ω, perpOf ((m.tbl i).v N) (m.vhat j N ω)⟫_ℝ|
          ≤ |⟪m.vhat i N ω, delocDir ((m.tbl i).v N) ((m.tbl j).X N ω)⟫_ℝ| :=
        abs_inner_perpOf_le ((m.tbl i).hv N) _ hω (m.mem_topSpace_vhat j N ω)
          (m.norm_vhat j N ω) _
      have h2 : |⟪m.vhat i N ω, delocDir ((m.tbl i).v N) ((m.tbl j).X N ω)⟫_ℝ|
          ≤ Real.sqrt (overlap ((m.tbl i).X N ω)
            (delocDir ((m.tbl i).v N) ((m.tbl j).X N ω))) := by
        rw [← Real.sqrt_sq_eq_abs]
        exact Real.sqrt_le_sqrt (overlap_ge_inner_sq _ _ (m.mem_topSpace_vhat i N ω)
          (m.norm_vhat i N ω))
      linarith
    · intro ε hε
      have hK := (law i).delocUniform (ε ^ 2) (by positivity)
      refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hK
        (fun _ => zero_le) (fun N => ?_)
      have hset : {ω : Ω N | ε ≤ |Real.sqrt (overlap ((m.tbl i).X N ω)
            (delocDir ((m.tbl i).v N) ((m.tbl j).X N ω))) - 0|}
          = {ω : Ω N | ε ^ 2 ≤ overlap ((m.tbl i).X N ω)
            (delocDir ((m.tbl i).v N) ((m.tbl j).X N ω))} := by
        ext ω
        simp only [sub_zero, Set.mem_ofPred_eq, abs_of_nonneg (Real.sqrt_nonneg _)]
        exact Real.le_sqrt hε.le (overlap_nonneg _ _)
      rw [hset]
      exact m.measure_deloc_le hI hij N (by positivity)
  have hsum := hA.add hB
  rw [add_zero] at hsum
  refine hsum.congr fun N => ?_
  filter_upwards with ω
  simp only [perpOf, inner_sub_right, real_inner_smul_right]
  rw [real_inner_comm (m.vhat j N ω) ((m.tbl i).v N)]
  ring

/-- The Gram matrix of `Ṽ` converges entrywise to `A_β`, signed: `⟪v̂_i, v̂_j⟫ → (A_β)_{ij}`.
The off-diagonal is `lem_delocalization`. The diagonal is exact, not a limit: `v̂_i` is a unit
vector because `SpikedModel.hd` gives `0 < d N`, so `(Ṽ Ṽᵀ)_{ii} = 1`, and
`(A_β)_{ii} = β_i² + (1 - β_i²) = 1`. The squared form cannot prove `thm_svd_stack_general`
(`notes/archive/audit_scope_2026-08-29.md` section 6.1). -/
theorem gram (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) (i j : Fin M) :
    TendstoInProb μ (fun N ω => ((m.Vt N ω) * (m.Vt N ω)ᵀ) i j)
      (Abeta (fun k => beta (m.tbl k).θ (c k)) i j) := by
  have hentry : ∀ (N : ℕ) (ω : Ω N),
      ⟪m.vhat i N ω, m.vhat j N ω⟫_ℝ = ((m.Vt N ω) * (m.Vt N ω)ᵀ) i j := by
    intro N ω
    rw [Matrix.mul_apply, real_inner_eq_dotProduct]
    rfl
  rcases eq_or_ne i j with rfl | hij
  · have hd : Abeta (fun k => beta (m.tbl k).θ (c k)) i i = 1 := by
      simp only [Abeta, Matrix.add_apply, Matrix.vecMulVec_apply, Matrix.diagonal_apply_eq]
      ring
    rw [hd]
    refine (TendstoInProb.const μ 1).congr fun N => ?_
    filter_upwards with ω
    rw [← hentry N ω, real_inner_self_eq_norm_sq, m.norm_vhat i N ω, one_pow]
  · have hd : Abeta (fun k => beta (m.tbl k).θ (c k)) i j
        = beta (m.tbl i).θ (c i) * beta (m.tbl j).θ (c j) := by
      simp only [Abeta, Matrix.add_apply, Matrix.vecMulVec_apply,
        Matrix.diagonal_apply_ne _ hij, add_zero]
    rw [hd]
    refine (m.lem_delocalization c law hI hij).congr fun N => ?_
    filter_upwards with ω
    exact hentry N ω

/-! ### The good event of the assembly

`G_N := Ṽ Ṽᵀ` is the `M × M` Gram matrix, `(G_N)_{ij} = ⟪v̂_i, v̂_j⟫`. The assembly of
`thm_svd_stack_general` is deterministic on the event where `∑_i P_i = Ṽᵀ Ṽ` and `G_N` has a
simple and positive top eigenvalue. That event has probability tending to one. -/

/-- The good event. -/
def goodEvent (m : MultiTableModel μ M n d) (N : ℕ) (ω : Ω N) : Prop :=
  m.svdstackGram N ω = (m.Vt N ω)ᵀ * m.Vt N ω ∧
    TopSimple (m.Vt N ω * (m.Vt N ω)ᵀ) (isHermitian_mul_transpose_self (m.Vt N ω)) ∧
    0 < lamMax (m.Vt N ω * (m.Vt N ω)ᵀ) (isHermitian_mul_transpose_self (m.Vt N ω))

/-- `gram` with the limit written through `β` instead of `beta (m.tbl k).θ (c k)`. -/
theorem gramEntries (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) (i j : Fin M) :
    TendstoInProb μ (fun N ω => (m.Vt N ω * (m.Vt N ω)ᵀ) i j) (Abeta β i j) := by
  have hβ : β = fun k => beta (m.tbl k).θ (c k) := funext hβdef
  have h := m.gram c law hI i j
  rw [← hβ] at h
  exact h

/-- The complement of the good event has vanishing probability. Route: `gramEntries` gives
the entrywise limit `G_N → A_β`; `topSimple_whp_of_tendsto` with the gap `abeta_gap` and
`lamMax_gt_half_whp_of_tendsto` with `λ_max(A_β) ≥ 1` (`one_le_lamMax_Abeta`) make the top
eigenvalue of `G_N` simple and positive with probability tending to one; `svdstackGram_eq` is
a.s. for each `N`. -/
theorem tendsto_measure_not_goodEvent (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ i j, i ≠ j ∧ 0 < β i ∧ 0 < β j)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) :
    Tendsto (fun N => μ N {ω | ¬ m.goodEvent N ω}) atTop (𝓝 0) := by
  obtain ⟨i, j, hij, hi, hj⟩ := hthr
  have h1 : 1 < Fintype.card (Fin M) := one_lt_card_fin_of_ne hij
  have hβ01 : ∀ k, 0 ≤ β k ∧ β k < 1 := fun k => by
    rw [hβdef k]
    exact beta_mem_Ico (hc k)
  have hA : (Abeta β).IsHermitian := isHermitian_Abeta β
  have hγ : 0 < lamMax (Abeta β) hA - hA.eigenvalues₀ ⟨1, h1⟩ :=
    lt_of_lt_of_le (mul_pos hi hj) (abeta_gap β hβ01 hij hi hj)
  have hApos : 0 < lamMax (Abeta β) hA := lt_of_lt_of_le zero_lt_one (one_le_lamMax_Abeta β i)
  have hsymm : ∀ (N : ℕ) (ω : Ω N), (m.Vt N ω * (m.Vt N ω)ᵀ).IsHermitian :=
    fun N ω => isHermitian_mul_transpose_self (m.Vt N ω)
  have hS := topSimple_whp_of_tendsto hA hsymm (m.gramEntries c β hβdef law hI) h1 hγ
  have hL := lamMax_gt_half_whp_of_tendsto hA hsymm (m.gramEntries c β hβdef law hI) hApos
  have hzero : ∀ N, μ N {ω | ¬ m.svdstackGram N ω = (m.Vt N ω)ᵀ * m.Vt N ω} = 0 :=
    fun N => ae_iff.mp (m.svdstackGram_eq c law N)
  have hsum : Tendsto (fun N => μ N {ω | ¬ TopSimple (m.Vt N ω * (m.Vt N ω)ᵀ) (hsymm N ω)} +
      μ N {ω | lamMax (m.Vt N ω * (m.Vt N ω)ᵀ) (hsymm N ω) ≤ lamMax (Abeta β) hA / 2})
      atTop (𝓝 0) := by
    simpa using hS.add hL
  have hsub : ∀ N, {ω | ¬ m.goodEvent N ω} ⊆
      {ω | ¬ m.svdstackGram N ω = (m.Vt N ω)ᵀ * m.Vt N ω} ∪
        ({ω | ¬ TopSimple (m.Vt N ω * (m.Vt N ω)ᵀ) (hsymm N ω)} ∪
          {ω | lamMax (m.Vt N ω * (m.Vt N ω)ᵀ) (hsymm N ω) ≤ lamMax (Abeta β) hA / 2}) := by
    intro N ω hω
    have hω' : ¬ m.goodEvent N ω := hω
    by_cases ha : m.svdstackGram N ω = (m.Vt N ω)ᵀ * m.Vt N ω
    · right
      by_cases hb : TopSimple (m.Vt N ω * (m.Vt N ω)ᵀ) (hsymm N ω)
      · right
        change lamMax (m.Vt N ω * (m.Vt N ω)ᵀ) (hsymm N ω) ≤ lamMax (Abeta β) hA / 2
        by_contra hc'
        rw [not_le] at hc'
        exact hω' ⟨ha, hb, by linarith⟩
      · left
        exact hb
    · left
      exact ha
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hsum
    (fun _ => zero_le) (fun N => ?_)
  calc μ N {ω | ¬ m.goodEvent N ω}
      ≤ μ N ({ω | ¬ m.svdstackGram N ω = (m.Vt N ω)ᵀ * m.Vt N ω} ∪
        ({ω | ¬ TopSimple (m.Vt N ω * (m.Vt N ω)ᵀ) (hsymm N ω)} ∪
          {ω | lamMax (m.Vt N ω * (m.Vt N ω)ᵀ) (hsymm N ω) ≤ lamMax (Abeta β) hA / 2})) :=
        measure_mono (hsub N)
    _ ≤ μ N {ω | ¬ m.svdstackGram N ω = (m.Vt N ω)ᵀ * m.Vt N ω} +
        μ N ({ω | ¬ TopSimple (m.Vt N ω * (m.Vt N ω)ᵀ) (hsymm N ω)} ∪
          {ω | lamMax (m.Vt N ω * (m.Vt N ω)ᵀ) (hsymm N ω) ≤ lamMax (Abeta β) hA / 2}) :=
        measure_union_le _ _
    _ ≤ μ N {ω | ¬ m.svdstackGram N ω = (m.Vt N ω)ᵀ * m.Vt N ω} +
        (μ N {ω | ¬ TopSimple (m.Vt N ω * (m.Vt N ω)ᵀ) (hsymm N ω)} +
          μ N {ω | lamMax (m.Vt N ω * (m.Vt N ω)ᵀ) (hsymm N ω) ≤ lamMax (Abeta β) hA / 2}) :=
        add_le_add le_rfl (measure_union_le _ _)
    _ = μ N {ω | ¬ TopSimple (m.Vt N ω * (m.Vt N ω)ᵀ) (hsymm N ω)} +
          μ N {ω | lamMax (m.Vt N ω * (m.Vt N ω)ᵀ) (hsymm N ω) ≤ lamMax (Abeta β) hA / 2} := by
        rw [hzero N, zero_add]

/-- On the good event the top eigenvalue of `∑_i P_i` is simple: `∑_i P_i = Ṽᵀ Ṽ`, whose
top eigenvalue is that of `G_N = Ṽ Ṽᵀ` (`lamMax_mul_transpose_self_eq`), and simplicity
transfers by `topSimple_transpose_mul_iff`. -/
private theorem topSimple_svdstackGram_of_goodEvent (m : MultiTableModel μ M n d) (N : ℕ)
    (ω : Ω N) (h : m.goodEvent N ω) :
    TopSimple (m.svdstackGram N ω) (m.isHermitian_svdstackGram N ω) := by
  obtain ⟨heq, hsimple, hpos⟩ := h
  rw [topSimple_congr_iff heq _ (isHermitian_transpose_mul_self (m.Vt N ω))]
  have hlam := lamMax_mul_transpose_self_eq (m.Vt N ω) hpos
  exact (topSimple_transpose_mul_iff (m.Vt N ω) _ _ hlam hpos).mpr hsimple

/-- The top eigenvalue of `∑_i P_i` is simple with probability tending to one, under the
hypotheses of `thm_svd_stack_general`. This is what turns the projector form `svdstackPerf`
into the paper's `|⟨v̂_svdstack, v⟩|²` and lets `svdstackPerf_eq_closed` apply. Gap 6 of
`notes/archive/audit_scope_2026-08-29.md` section 3.2.

The statement is not almost sure for each `N`. For `M = 2`, `d_N = 2` simplicity fails on
`{⟪v̂_1, v̂_2⟫ = 0}`, a nonempty event; to call it null one needs a density for the joint law
of the top eigenvectors `(v̂_1, ..., v̂_M)`, which no hypothesis gives (`SingleTableLaw` is a
list of limits and one a.s. fact per table, and `∑_i P_i` is not a polynomial in the noise).
The final theorem is a convergence in probability, so it needs only this form
(user decision 5 of `notes/archive/thm_svd_stack_general.md`, 2026-08-29). -/
theorem topSimple_svdstackGram (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ i j, i ≠ j ∧ 0 < β i ∧ 0 < β j)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) :
    Tendsto (fun N => μ N {ω | ¬ TopSimple (m.svdstackGram N ω)
      (m.isHermitian_svdstackGram N ω)}) atTop (𝓝 0) := by
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds
    (m.tendsto_measure_not_goodEvent c β hc hβdef hthr law hI) (fun _ => zero_le)
    (fun N => measure_mono fun ω hω => ?_)
  have hω' : ¬ TopSimple (m.svdstackGram N ω) (m.isHermitian_svdstackGram N ω) := hω
  change ¬ m.goodEvent N ω
  intro hgood
  exact hω' (m.topSimple_svdstackGram_of_goodEvent N ω hgood)

/-! ### The rows of `Ṽ v` -/

-- One table at least. `Fin M` must be nonempty: the shared `v` is read off table `0`.
variable [NeZero M]

/-- Row `i` of `Ṽ v` converges to `β_i`, signed: `⟪v̂_i, v⟫ → β_i`. This is
`prop:single_table` part 1 for table `i`, plus the sign convention of `vhat`, which forces
`⟪v̂_i, v⟫ ≥ 0` and so picks the nonnegative square root. -/
theorem align (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (i : Fin M) :
    TendstoInProb μ (fun N ω => (m.Vt N ω).mulVec (WithLp.ofLp ((m.tbl 0).v N)) i)
      (beta (m.tbl i).θ (c i)) := by
  refine (m.align_inner c law i).congr fun N => ?_
  filter_upwards with ω
  rw [m.hv 0 i N, real_inner_eq_dotProduct]
  rfl

end SVDStack

end MultiTableModel

end StackedSVD
