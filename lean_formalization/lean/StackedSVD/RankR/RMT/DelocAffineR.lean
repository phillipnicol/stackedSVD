/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.DelocR
import StackedSVD.RMT.Symmetry
import StackedSVD.RankR.RMT.SimplicityAffineR
import StackedSVD.RankR.RMT.TableStack
import StackedSVD.LinAlg.SpecIdxMeas

/-!
# Task C6: `delocUniform` of `TableLawR` at general `rk`

Task C6 of `notes/archive/rankr_plan_C.md` section 2. The file proves the `delocUniform` field of
`SpikedModelR.TableLawR` (`RankR/General.lean:163`) for Gaussian noise, in both regimes, with no
hypothesis on `θ`. The rank-1 mirror is `SpikedModel.delocUniform_of_gaussian`
(`RMT/Symmetry.lean:603`), and the route is the rank-`rk` port of `RMT/Symmetry.lean` sections 7
to 10, with the top projector replaced by the projector at one eigenvalue index and the single
protected direction replaced by the `rk` columns of `V`.

## The argument

Write `X = U Θ Vᵀ + t Z` with `Z` a canonical Gaussian matrix and `t = d^{-1/2}`.

1. **Right rotation moves the test direction** (deterministic). For `Oᵀ O = 1` the Gram of
   `X O` is `Oᵀ (Xᵀ X) O`, and `DelocR.normSq_specProjIdx_mul_right` moves the index projector,
   so `overlapIdx (X O) k w = overlapIdx X k (O w)` (`overlapIdx_mul_right`).
2. **Right rotation preserves the law** (`Symmetry.measurePreserving_mul_right`) and is a
   measurable embedding (`DelocR.measurableEmbedding_mul_right`), so the change of variables
   asks nothing of the integrand.
3. **The signal is fixed.** `U Θ Vᵀ O = U Θ Vᵀ` as soon as `Oᵀ` fixes every column of `V`
   (`mul_right_of_transpose_mulVec_col`).
4. **Exchangeability** (`lintegral_overlapIdx_eq_of_orth`). For two unit vectors `w` and `w'`
   orthogonal to every column of `V`, the Householder reflection with axis `w - w'` fixes every
   such column, because both `w` and `w'` are orthogonal to it. Steps 1 to 3 then give equal
   integrals. This is risk 6 of the plan, on the primary route and not on the fallback.
5. **Bessel** (`sum_overlapIdx_le_one`). On the almost sure event `SimpleSpec _ (k+1)` of
   `SimplicityAffineR.simpleSpec_ae_affine` the projector at the index `k` is rank one, so the
   sum over an orthonormal family is at most `1` (`DelocR.sum_normSq_specProjIdx_le_one`).
6. **Count.** The family is an orthonormal basis of the orthogonal complement of the column span
   of `V`. The columns are orthonormal (`hV`), so that span has dimension `rk` and the
   complement has dimension `p - rk`; the bound is `1/(p - rk)`.
7. **Markov** turns the mean bound into a tail bound, uniformly over `w ∈ orthUnitR N`, and
   `d N → ∞` finishes.

Nothing here needs `Regime` beyond `d N → ∞`, and nothing needs `hθanti`, `hθnn` or a threshold
on `θ`; section 6 of the plan lists `delocUniform` among the results that stay true at a tie.

## The strict rank bound

Minor finding 5 of `notes/archive/audit_rankr_plan_C_2026-09-02.md`: the mean bound is `1/(d - rk)`,
so `delocUniform_of_gaussian` needs `rk < d N` and not the `rk ≤ d N` of `SpikedModelR.rk_le_d`. The
target is a limit, so `Filter.Eventually` is enough, and `hd : Tendsto d atTop atTop` supplies `rk +
1 ≤ d N` for large `N`. `rk_le_n` and `rk_le_d` (`RankR/RMT/TableStack.lean:111,116`) supply the
index bound `k < min (n N) (d N)` at every `N`.

## Measurability

The mean bound of section 1 needs none: the change of variables goes through
`MeasurePreserving.lintegral_comp_emb` and the Bessel step through `DelocR.sum_lintegral_le`,
and neither asks anything of the integrand. Markov does need it, and it comes from
`measurable_overlapIdx` (`LinAlg/SpecIdxMeas.lean:411`) composed with the affine map. The
rank-1 file keeps that composition `private` (`RMT/Symmetry.lean:368,376`), so this file
restates both helpers privately, as attack 4 of the audit asks of C7.

No `sorry`, no `axiom`, no edit to any existing file.
-/

open MeasureTheory Filter Topology
open scoped InnerProductSpace Matrix ENNReal

namespace StackedSVD

/-! ### 0. One private helper restated from `RMT/Symmetry.lean`

`measurable_affine_matrix` and `overlapIdx_nonneg` used to stand here as private copies. The
second cleanup pass (2026-09-02) moved one public copy of each into `LinAlg/SpecIdxMeas.lean`
and `LinAlg/SpecIdx.lean`. -/

/-- Index form of the `private` `measurable_overlap_affine` (`RMT/Symmetry.lean:376`), from
`measurable_overlapIdx` (`LinAlg/SpecIdxMeas.lean`) and `measurable_affine_matrix`
(`LinAlg/SpecIdxMeas.lean`). -/
private theorem measurable_overlapIdx_affine {q p : ℕ} (A : Matrix (Fin q) (Fin p) ℝ) (t : ℝ)
    (kk : ℕ) (w : EuclideanSpace ℝ (Fin p)) :
    Measurable fun Z : Matrix (Fin q) (Fin p) ℝ =>
      ENNReal.ofReal (overlapIdx (A + t • Z) kk w) :=
  ENNReal.measurable_ofReal.comp
    ((measurable_overlapIdx kk w).comp (measurable_affine_matrix A t))

/-! ### 1. The three deterministic steps -/

/-- **Right multiplication by an orthogonal matrix moves the test direction.** Index form of
`Symmetry.overlap_mul_right`; at `kk = 0` the two agree. -/
theorem overlapIdx_mul_right {q p : ℕ} (X : Matrix (Fin q) (Fin p) ℝ)
    {O : Matrix (Fin p) (Fin p) ℝ} (hO : Oᵀ * O = 1) (kk : ℕ) (w : EuclideanSpace ℝ (Fin p)) :
    overlapIdx (X * O) kk w = overlapIdx X kk (rotIso O hO w) := by
  rw [overlapIdx, overlapIdx, normSq_specProjIdx_mul_right X hO kk w]

/-- **The rank-`rk` signal is fixed by a right rotation that fixes every column of `V`.** The
rank-`rk` twin of `Symmetry.vecMulVec_mul_of_transpose_mulVec_eq`, which is the case `rk = 1`
written with `vecMulVec`. -/
theorem mul_right_of_transpose_mulVec_col {q p rk : ℕ} (U : Matrix (Fin q) (Fin rk) ℝ)
    (θ : Fin rk → ℝ) (V : Matrix (Fin p) (Fin rk) ℝ) {O : Matrix (Fin p) (Fin p) ℝ}
    (hO : ∀ a : Fin rk, Oᵀ *ᵥ (fun l => V l a) = fun l => V l a) :
    U * Matrix.diagonal θ * Vᵀ * O = U * Matrix.diagonal θ * Vᵀ := by
  have hcol : Oᵀ * V = V := by
    ext l a
    rw [Matrix.mul_apply, ← Matrix.mulVec_apply_eq_sum]
    exact congrFun (hO a) l
  have hVO : Vᵀ * O = Vᵀ := by
    calc Vᵀ * O = (Oᵀ * V)ᵀ := by rw [Matrix.transpose_mul, Matrix.transpose_transpose]
      _ = Vᵀ := by rw [hcol]
  rw [Matrix.mul_assoc, hVO]

/-! ### 2. Bessel at one index -/

/-- **Bessel at one index, in overlap form.** On the event that the eigenvalue at the sorted
index `kk` is simple, the projector there is rank one, so the overlaps with an orthonormal
family sum to at most `1`. Index form of `Symmetry.sum_overlap_le_one`. -/
theorem sum_overlapIdx_le_one {q p : ℕ} {ι : Type*} [Fintype ι] {X : Matrix (Fin q) (Fin p) ℝ}
    {kk : ℕ} (hsimple : SimpleSpec (Xᵀ * X) (isHermitian_transpose_mul_self X) (kk + 1))
    (hk : kk < p) {e : ι → EuclideanSpace ℝ (Fin p)} (he : Orthonormal ℝ e) :
    ∑ i, overlapIdx X kk (e i) ≤ 1 :=
  sum_normSq_specProjIdx_le_one hsimple hk he

/-! ### 3. Exchangeability of the test direction -/

/-- **The law of `overlapIdx X kk w` is the same for every unit `w` orthogonal to every column
of `V`.** The rank-`rk` twin of `Symmetry.lintegral_overlap_eq_of_orth`. The Householder
reflection with axis `w - w'` fixes every column of `V`, because both `w` and `w'` are
orthogonal to it, so the signal invariance hypothesis `hA` applies to it. -/
theorem lintegral_overlapIdx_eq_of_orth {q p rk : ℕ} (A : Matrix (Fin q) (Fin p) ℝ) (t : ℝ)
    (V : Matrix (Fin p) (Fin rk) ℝ) (kk : ℕ)
    (hA : ∀ O : Matrix (Fin p) (Fin p) ℝ,
      (∀ a : Fin rk, Oᵀ *ᵥ (fun l => V l a) = fun l => V l a) → A * O = A)
    {w w' : EuclideanSpace ℝ (Fin p)} (hw : ‖w‖ = 1) (hw' : ‖w'‖ = 1)
    (hwV : ∀ a : Fin rk,
      ⟪w, (WithLp.toLp 2 fun l => V l a : EuclideanSpace ℝ (Fin p))⟫_ℝ = 0)
    (hw'V : ∀ a : Fin rk,
      ⟪w', (WithLp.toLp 2 fun l => V l a : EuclideanSpace ℝ (Fin p))⟫_ℝ = 0) :
    ∫⁻ Z, ENNReal.ofReal (overlapIdx (A + t • Z) kk w) ∂(gaussianMatrix q p)
      = ∫⁻ Z, ENNReal.ofReal (overlapIdx (A + t • Z) kk w') ∂(gaussianMatrix q p) := by
  by_cases hww : w = w'
  · rw [hww]
  have hdot : ∀ x y : EuclideanSpace ℝ (Fin p),
      WithLp.ofLp x ⬝ᵥ WithLp.ofLp y = ⟪x, y⟫_ℝ := fun x y =>
    (inner_euclidean_eq_dotProduct x y).symm
  have hww1 : WithLp.ofLp w ⬝ᵥ WithLp.ofLp w = 1 := by
    rw [hdot, real_inner_self_eq_norm_sq, hw]
    norm_num
  have hww'1 : WithLp.ofLp w' ⬝ᵥ WithLp.ofLp w' = 1 := by
    rw [hdot, real_inner_self_eq_norm_sq, hw']
    norm_num
  have hne : WithLp.ofLp w ≠ WithLp.ofLp w' := fun h => hww (by
    have := congrArg (WithLp.toLp 2) h
    simpa using this)
  set ax : Fin p → ℝ := WithLp.ofLp w - WithLp.ofLp w' with hax
  set O : Matrix (Fin p) (Fin p) ℝ := householder ax with hOdef
  have hax0 : ax ⬝ᵥ ax ≠ 0 := fun h => (sub_ne_zero.mpr hne) (dotProduct_self_eq_zero.mp h)
  have hO : Oᵀ * O = 1 := householder_orth hax0
  -- the axis is orthogonal to every column of `V`, so the reflection fixes every column
  have hfix : ∀ a : Fin rk, Oᵀ *ᵥ (fun l => V l a) = fun l => V l a := by
    intro a
    have h1 : WithLp.ofLp w ⬝ᵥ (fun l => V l a) = 0 := by
      have h := hdot w (WithLp.toLp 2 fun l => V l a)
      rw [hwV a] at h
      simpa using h
    have h2 : WithLp.ofLp w' ⬝ᵥ (fun l => V l a) = 0 := by
      have h := hdot w' (WithLp.toLp 2 fun l => V l a)
      rw [hw'V a] at h
      simpa using h
    have hac : ax ⬝ᵥ (fun l => V l a) = 0 := by
      rw [hax, sub_dotProduct, h1, h2, sub_zero]
    rw [hOdef, householder_transpose]
    exact householder_apply_of_orth hac
  have hOw : rotIso O hO w = w' := by
    rw [rotIso_apply]
    exact congrArg (WithLp.toLp 2) (householder_sub_apply hww1 hww'1 hne)
  -- the change of variables
  have key : ∀ Z : Matrix (Fin q) (Fin p) ℝ,
      overlapIdx (A + t • (Z * O)) kk w = overlapIdx (A + t • Z) kk w' := by
    intro Z
    have hmul : (A + t • Z) * O = A + t • (Z * O) := by
      rw [Matrix.add_mul, hA O hfix, Matrix.smul_mul]
    rw [← hmul, overlapIdx_mul_right _ hO, hOw]
  calc ∫⁻ Z, ENNReal.ofReal (overlapIdx (A + t • Z) kk w) ∂(gaussianMatrix q p)
      = ∫⁻ Z, ENNReal.ofReal (overlapIdx (A + t • (Z * O)) kk w) ∂(gaussianMatrix q p) :=
        ((measurePreserving_mul_right hO q).lintegral_comp_emb
          (measurableEmbedding_mul_right hO q) _).symm
    _ = _ := by simp only [key]

/-! ### 4. The mean bound `1/(p - rk)` -/

/-- **Item C6.1.** For the affine family `U Θ Vᵀ + t Z` with `V` of orthonormal columns, for a
sorted eigenvalue index `kk < min q p` and for a unit direction `w` orthogonal to every column
of `V`, the mean overlap at that index is at most `1/(p - rk)`.

The rank-1 mirror is `Symmetry.lintegral_overlap_le` (`RMT/Symmetry.lean:467`), whose bound is
`1/(p - 1)`: one protected direction there, `rk` of them here. -/
theorem lintegral_overlapIdx_affine_le {q p : ℕ} (hq : 0 < q) {rk : ℕ} (hrk : rk < p)
    (U : Matrix (Fin q) (Fin rk) ℝ) (θ : Fin rk → ℝ) (V : Matrix (Fin p) (Fin rk) ℝ)
    (hV : Vᵀ * V = 1) {t : ℝ} (ht : t ≠ 0) {kk : ℕ} (hk : kk < min q p)
    {w : EuclideanSpace ℝ (Fin p)} (hw : ‖w‖ = 1)
    (hwV : ∀ a : Fin rk,
      ⟪w, (WithLp.toLp 2 fun l => V l a : EuclideanSpace ℝ (Fin p))⟫_ℝ = 0) :
    ∫⁻ Z, ENNReal.ofReal (overlapIdx (U * Matrix.diagonal θ * Vᵀ + t • Z) kk w)
        ∂(gaussianMatrix q p) ≤ ENNReal.ofReal (1 / ((p : ℝ) - rk)) := by
  have hp : 0 < p := lt_of_le_of_lt (Nat.zero_le _) hrk
  set A : Matrix (Fin q) (Fin p) ℝ := U * Matrix.diagonal θ * Vᵀ with hAdef
  have hAinv : ∀ O : Matrix (Fin p) (Fin p) ℝ,
      (∀ a : Fin rk, Oᵀ *ᵥ (fun l => V l a) = fun l => V l a) → A * O = A := by
    intro O hO
    rw [hAdef]
    exact mul_right_of_transpose_mulVec_col U θ V hO
  -- the columns of `V`, as an orthonormal family
  set v : Fin rk → EuclideanSpace ℝ (Fin p) :=
    fun a => (WithLp.toLp 2 fun l => V l a : EuclideanSpace ℝ (Fin p)) with hvdef
  have hcolinner : ∀ a b : Fin rk, ⟪v a, v b⟫_ℝ = (Vᵀ * V) a b := by
    intro a b
    rw [inner_euclidean_eq_dotProduct, Matrix.mul_apply]
    simp [hvdef, dotProduct, Matrix.transpose_apply]
  have hvon : Orthonormal ℝ v := by
    rw [orthonormal_iff_ite]
    intro a b
    rw [hcolinner a b, hV]
    simp [Matrix.one_apply]
  have hwv : ∀ a : Fin rk, ⟪w, v a⟫_ℝ = 0 := hwV
  -- an orthonormal basis of the orthogonal complement of the column span
  set W : Submodule ℝ (EuclideanSpace ℝ (Fin p)) := Submodule.span ℝ (Set.range v) with hWdef
  have hWrank : Module.finrank ℝ W = rk := by
    rw [hWdef, finrank_span_eq_card hvon.linearIndependent, Fintype.card_fin]
  have hOrank : Module.finrank ℝ (Wᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin p))) = p - rk := by
    have h2 := Submodule.finrank_add_finrank_orthogonal (𝕜 := ℝ) (K := W)
    rw [hWrank, finrank_euclideanSpace_fin] at h2
    omega
  set b := stdOrthonormalBasis ℝ (Wᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin p))) with hb
  set e : Fin (Module.finrank ℝ (Wᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin p)))) →
      EuclideanSpace ℝ (Fin p) := fun j => (b j : EuclideanSpace ℝ (Fin p)) with he
  have hon : Orthonormal ℝ e :=
    ((Wᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin p))).subtypeₗᵢ.orthonormal_comp_iff
      (v := fun j => b j)).mpr b.orthonormal
  have henorm : ∀ j, ‖e j‖ = 1 := fun j => hon.1 j
  have hev : ∀ j, ∀ a : Fin rk, ⟪e j, v a⟫_ℝ = 0 := by
    intro j a
    have hmem : (b j : EuclideanSpace ℝ (Fin p)) ∈ (Wᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin p))) :=
      (b j).2
    have hva : v a ∈ W := Submodule.subset_span ⟨a, rfl⟩
    have hzero := (Submodule.mem_orthogonal W _).mp hmem (v a) hva
    rw [real_inner_comm]
    exact hzero
  -- every admissible direction has the same mean overlap
  set I : ℝ≥0∞ := ∫⁻ Z, ENNReal.ofReal (overlapIdx (A + t • Z) kk w) ∂(gaussianMatrix q p)
    with hI
  have hsame : ∀ j, ∫⁻ Z, ENNReal.ofReal (overlapIdx (A + t • Z) kk (e j))
      ∂(gaussianMatrix q p) = I := fun j =>
    (lintegral_overlapIdx_eq_of_orth A t V kk hAinv hw (henorm j) hwv (hev j)).symm
  -- the tie set is a null set: on its complement the projector at the index `kk` is rank one
  have hae : ∀ᵐ Z ∂(gaussianMatrix q p),
      ∑ j, ENNReal.ofReal (overlapIdx (A + t • Z) kk (e j)) ≤ 1 := by
    filter_upwards [simpleSpec_ae_affine hq hp A ht (kk + 1) (by omega)] with Z hZ
    have hbes : ∑ j, overlapIdx (A + t • Z) kk (e j) ≤ 1 :=
      sum_overlapIdx_le_one hZ (by omega) hon
    calc ∑ j, ENNReal.ofReal (overlapIdx (A + t • Z) kk (e j))
        = ENNReal.ofReal (∑ j, overlapIdx (A + t • Z) kk (e j)) :=
          (ENNReal.ofReal_sum_of_nonneg fun j _ => overlapIdx_nonneg _ _ _).symm
      _ ≤ 1 := ENNReal.ofReal_le_one.mpr hbes
  have hsum : ∑ j, ∫⁻ Z, ENNReal.ofReal (overlapIdx (A + t • Z) kk (e j))
      ∂(gaussianMatrix q p) ≤ 1 := by
    calc ∑ j, ∫⁻ Z, ENNReal.ofReal (overlapIdx (A + t • Z) kk (e j)) ∂(gaussianMatrix q p)
        ≤ ∫⁻ Z, ∑ j, ENNReal.ofReal (overlapIdx (A + t • Z) kk (e j))
            ∂(gaussianMatrix q p) := sum_lintegral_le _ _
      _ ≤ ∫⁻ _, (1 : ℝ≥0∞) ∂(gaussianMatrix q p) := lintegral_mono_ae hae
      _ = 1 := by simp
  -- conclude
  simp only [hsame, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at hsum
  have hcard : (Module.finrank ℝ (Wᗮ : Submodule ℝ (EuclideanSpace ℝ (Fin p))) : ℝ≥0∞)
      = ENNReal.ofReal ((p : ℝ) - rk) := by
    rw [hOrank]
    rw [show ((p - rk : ℕ) : ℝ≥0∞) = ENNReal.ofReal ((p - rk : ℕ) : ℝ) by
      rw [ENNReal.ofReal_natCast]]
    congr 1
    have hle : rk ≤ p := hrk.le
    push_cast [Nat.cast_sub hle]
    ring
  rw [hcard] at hsum
  have hpos : (0 : ℝ) < (p : ℝ) - rk := by
    have : (rk : ℝ) < (p : ℝ) := by exact_mod_cast hrk
    linarith
  have hne0 : ENNReal.ofReal ((p : ℝ) - rk) ≠ 0 := by
    simp only [ne_eq, ENNReal.ofReal_eq_zero, not_le]
    exact hpos
  have hnetop : ENNReal.ofReal ((p : ℝ) - rk) ≠ ⊤ := ENNReal.ofReal_ne_top
  have hdiv : I ≤ 1 / ENNReal.ofReal ((p : ℝ) - rk) := by
    rw [ENNReal.le_div_iff_mul_le (Or.inl hne0) (Or.inl hnetop), mul_comm]
    exact hsum
  refine hdiv.trans (le_of_eq ?_)
  rw [ENNReal.ofReal_div_of_pos hpos, ENNReal.ofReal_one]

/-! ### 5. Transfer to the model, Markov, and `delocUniform` -/

namespace SpikedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}
  {rk : ℕ}

/-- The mean overlap bound for the model, at one `N`. The rank-1 mirror is
`SpikedModel.lintegral_overlap_model_le` (`RMT/Symmetry.lean:546`). -/
theorem lintegral_overlapIdx_model_le (m : SpikedModelR μ n d rk) (hG : m.GaussianNoise)
    (N : ℕ) (hrk : rk < d N) (k : Fin rk) {w : EuclideanSpace ℝ (Fin (d N))}
    (hw : w ∈ m.orthUnitR N) :
    ∫⁻ ω, ENNReal.ofReal (overlapIdx (m.X N ω) (k : ℕ) w) ∂(μ N)
      ≤ ENNReal.ofReal (1 / ((d N : ℝ) - rk)) := by
  set A : Matrix (Fin (n N)) (Fin (d N)) ℝ := m.U N * Matrix.diagonal m.θ * (m.V N)ᵀ with hAdef
  set t : ℝ := (Real.sqrt (d N))⁻¹ with htdef
  have hdpos : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast m.hd N
  have ht : t ≠ 0 := inv_ne_zero (Real.sqrt_pos.mpr hdpos).ne'
  have hk : (k : ℕ) < min (n N) (d N) :=
    lt_min (lt_of_lt_of_le k.isLt (m.rk_le_n N)) (lt_of_lt_of_le k.isLt (m.rk_le_d N))
  have hX : ∀ ω, m.X N ω = A + t • m.Z N ω := fun ω => rfl
  have hfun : (fun Z : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
      ENNReal.ofReal (overlapIdx (A + t • Z) (k : ℕ) w)) ∘ (m.Z N) =
      fun ω => ENNReal.ofReal (overlapIdx (m.X N ω) (k : ℕ) w) := by
    funext ω
    rw [Function.comp_apply, hX ω]
  calc ∫⁻ ω, ENNReal.ofReal (overlapIdx (m.X N ω) (k : ℕ) w) ∂(μ N)
      = ∫⁻ Z, ENNReal.ofReal (overlapIdx (A + t • Z) (k : ℕ) w)
          ∂(gaussianMatrix (n N) (d N)) := by
        rw [← hfun]
        exact (hG N).lintegral_comp (measurable_overlapIdx_affine A t (k : ℕ) w).aemeasurable
    _ ≤ ENNReal.ofReal (1 / ((d N : ℝ) - rk)) :=
        lintegral_overlapIdx_affine_le (m.hn N) hrk (m.U N) m.θ (m.V N) (m.hV N) ht hk hw.1 hw.2

/-- Markov's inequality turns the mean bound into a tail bound. The rank-1 mirror is
`SpikedModel.measure_overlap_ge_le` (`RMT/Symmetry.lean:573`). -/
theorem measure_overlapIdx_ge_le (m : SpikedModelR μ n d rk) (hG : m.GaussianNoise) (N : ℕ)
    (hrk : rk < d N) (k : Fin rk) {ε : ℝ} (hε : 0 < ε)
    {w : EuclideanSpace ℝ (Fin (d N))} (hw : w ∈ m.orthUnitR N) :
    μ N {ω | ε ≤ overlapIdx (m.X N ω) (k : ℕ) w}
      ≤ ENNReal.ofReal (1 / (ε * ((d N : ℝ) - rk))) := by
  have hmeas : AEMeasurable
      (fun ω => ENNReal.ofReal (overlapIdx (m.X N ω) (k : ℕ) w)) (μ N) := by
    have h1 : Measurable fun ω => ENNReal.ofReal (overlapIdx (m.X N ω) (k : ℕ) w) := by
      have hcomp : (fun ω => ENNReal.ofReal (overlapIdx (m.X N ω) (k : ℕ) w)) =
          (fun Z : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
            ENNReal.ofReal (overlapIdx (m.U N * Matrix.diagonal m.θ * (m.V N)ᵀ
              + (Real.sqrt (d N))⁻¹ • Z) (k : ℕ) w)) ∘ (m.Z N) := rfl
      rw [hcomp]
      exact (measurable_overlapIdx_affine _ _ _ w).comp (m.hZ N)
    exact h1.aemeasurable
  have hset : {ω | ε ≤ overlapIdx (m.X N ω) (k : ℕ) w}
      = {ω | ENNReal.ofReal ε ≤ ENNReal.ofReal (overlapIdx (m.X N ω) (k : ℕ) w)} := by
    ext ω
    simp [ENNReal.ofReal_le_ofReal_iff (overlapIdx_nonneg _ _ _)]
  have hεne : ENNReal.ofReal ε ≠ 0 := by simp [ENNReal.ofReal_eq_zero, not_le, hε]
  have hmk := meas_ge_le_lintegral_div hmeas hεne ENNReal.ofReal_ne_top
  rw [hset]
  refine hmk.trans ?_
  have hbound := lintegral_overlapIdx_model_le m hG N hrk k hw
  refine (ENNReal.div_le_div_right hbound _).trans (le_of_eq ?_)
  rw [← ENNReal.ofReal_div_of_pos hε]
  congr 1
  field_simp

/-- **Item C6.2: the `delocUniform` field of `TableLawR`, for Gaussian noise.** Both regimes, no
hypothesis on `θ`, and uniform over `w ∈ orthUnitR N` with no extra work. The rank-1 mirror is
`SpikedModel.delocUniform_of_gaussian` (`RMT/Symmetry.lean:603`).

`hd` enters twice: it gives `rk < d N` for large `N` (the bound is `1/(d - rk)`, so `rk ≤ d N`
is not enough), and it drives the bound to `0`. -/
theorem delocUniform_of_gaussian (m : SpikedModelR μ n d rk) (hG : m.GaussianNoise)
    (hd : Tendsto d atTop atTop) (k : Fin rk) :
    ∀ ε > 0, Tendsto (fun N =>
      ⨆ w ∈ m.orthUnitR N, μ N {ω | ε ≤ overlapIdx (m.X N ω) (k : ℕ) w}) atTop (𝓝 0) := by
  intro ε hε
  have hle : ∀ᶠ N in atTop,
      (⨆ w ∈ m.orthUnitR N, μ N {ω | ε ≤ overlapIdx (m.X N ω) (k : ℕ) w})
        ≤ ENNReal.ofReal (1 / (ε * ((d N : ℝ) - rk))) := by
    filter_upwards [hd (Filter.eventually_ge_atTop (rk + 1))] with N hN
    have hN' : rk + 1 ≤ d N := hN
    exact iSup₂_le fun w hw => measure_overlapIdx_ge_le m hG N (by omega) k hε hw
  have hreal : Tendsto (fun N => 1 / (ε * ((d N : ℝ) - rk))) atTop (𝓝 0) := by
    have h1 : Tendsto (fun N => ε * ((d N : ℝ) - rk)) atTop atTop := by
      refine Filter.Tendsto.const_mul_atTop hε ?_
      exact (tendsto_natCast_atTop_atTop.comp hd).atTop_add tendsto_const_nhds
    exact Filter.Tendsto.congr (fun N => (one_div _).symm) h1.inv_tendsto_atTop
  have htend : Tendsto (fun N => ENNReal.ofReal (1 / (ε * ((d N : ℝ) - rk)))) atTop (𝓝 0) := by
    have h3 := (ENNReal.continuous_ofReal.tendsto 0).comp hreal
    rw [ENNReal.ofReal_zero] at h3
    exact h3
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds htend ?_ hle
  exact Filter.Eventually.of_forall fun _ => by simp

end SpikedModelR

end StackedSVD
