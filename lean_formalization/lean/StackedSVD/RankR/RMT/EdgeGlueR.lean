/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.EdgeGlueDetR
import StackedSVD.RankR.RMT.EdgeR
import StackedSVD.RankR.RMT.SymmetryR
import StackedSVD.RankR.RMT.SimplicityAffineR
import StackedSVD.LinAlg.SpecIdxMeas

/-!
# Task U7c, gap G1 glue: the probabilistic layer

Gap G1 of `notes/archive/rankr_plan_A.md` section 4. The overlap of a spike direction with the
top-`r` projector of the Gram matrix splits (`Frame.norm_sq_specProjTop_split`) into the part
carried by eigenvalues above `bulkEdge c + ε₁` and the part carried by the top-`r` eigenvalues at or
below `bulkEdge c + ε₁`. This file bounds the second part: for every spike `j` and every `δ > 0`
there is an `ε₀ > 0` such that every `ε₁ ∈ (0, ε₀]` gives that part at most `δ`, with probability
tending to `1`.

No supercriticality is assumed, so the statement holds for every spike of the mixed regime.

## Route

`EdgeGlueDetR.normSq_specProj_edge_le` does the linear algebra at one sample point. This file
produces the seven good events it consumes, each of probability tending to `1`, and takes a
union bound in the style of `RMT/R6.lean:512-757`.

| event | statement | source |
|---|---|---|
| G1 | `lamMax W₀ ≤ b + ε₀` and `lamMax W₀ ≤ b + 1/2` | `ResolventLimitsR.edge` |
| G2 | `SimpleSpec W₀ r` | `simpleSpec_ae_affine` at `A = 0` |
| G3 | `SimpleSpec (Xᵀ X) r` | `simpleSpec_ae_affine` at `A = signalPart` |
| G4 | the κ bound on the top-`r` projectors of `W₀` at the columns of `Q` | `DelocR` and Markov |
| G5 | the `r²` entries of `Qᵀ G₀(z₀)² Q` near `diag((λ + 1) m')` | `tendstoInProb_cform2_qmatR` |
| G6 | the `r²` entries of `Qᵀ G₀(z₁) Q` near `diag((λ + 1) m₁)` | `tendstoInProb_cform_qmatR` |
| G7 | the residual direction has small top-`r` overlap | `SymmetryR` and Markov |

The size conditions `r ≤ pp N` and `r < d N` hold for large `N` only, so every inclusion is an
eventual one.

## The truncation in G4

`DelocR.lintegral_prod_normSq_specProjIdx_le` bounds the mean of `‖P_a q‖²` by
`(∫⁻ ‖q‖²) / d`, and the column `q` of `Q` has an unbounded Gaussian part, whose second moment
this project does not compute anywhere. The κ event is therefore stated on the truncated column
`if ‖q‖² ≤ C then q else 0`, whose norm is bounded by construction. Event G5 supplies
`‖q‖² ≤ C` at the same sample point, so the truncation is invisible on the intersection of the
good events. `C` comes from `EdgeGlueDetR.le_qform2_of_psd`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

namespace EdgeGlueR

/-! ### 1. Small measurability and congruence helpers -/

section Helpers

/-- `SimpleSpec` moves along an equality of matrices. -/
theorem simpleSpec_congr {D : ℕ} {A A' : Matrix (Fin D) (Fin D) ℝ} (h : A = A')
    (hA : A.IsHermitian) (hA' : A'.IsHermitian) {rk : ℕ} (hs : SimpleSpec A hA rk) :
    SimpleSpec A' hA' rk := by
  subst h
  exact hs

/-- `‖Q y‖² = yᵀ QᵀQ y`. -/
theorem dotProduct_mulVec_self {D rr : ℕ} (Q : Matrix (Fin D) (Fin rr) ℝ) (y : Fin rr → ℝ) :
    (Q *ᵥ y) ⬝ᵥ (Q *ᵥ y) = y ⬝ᵥ ((Qᵀ * Q) *ᵥ y) := by
  rw [EdgeGlueDetR.dotProduct_mulVec_transpose Q (Q *ᵥ y) y, Matrix.mulVec_mulVec,
    dotProduct_comm]

/-- The determinant of a matrix whose entries are measurable is measurable. -/
theorem measurable_det_of_entries {X : Type*} [MeasurableSpace X] {n : ℕ}
    {M : X → Matrix (Fin n) (Fin n) ℝ} (hM : ∀ i j, Measurable fun x => M x i j) :
    Measurable fun x => (M x).det := by
  classical
  have h : (fun x => (M x).det) = fun x => ∑ σ : Equiv.Perm (Fin n),
      ((Equiv.Perm.sign σ : ℤ) : ℝ) * ∏ i, M x (σ i) i := by
    funext x
    exact Matrix.det_apply' (M x)
  rw [h]
  refine Finset.measurable_sum _ fun σ _ => ?_
  exact measurable_const.mul (Finset.measurable_prod _ fun i _ => hM (σ i) i)

/-- The set of matrices whose Gram matrix has a simple top-`rk` spectrum is measurable. -/
theorem measurableSet_simpleSpec_gram {n D : ℕ} (rk : ℕ) :
    MeasurableSet {Y : Matrix (Fin n) (Fin D) ℝ |
      SimpleSpec (Yᵀ * Y) (isHermitian_transpose_mul_self Y) rk} := by
  classical
  have hset : {Y : Matrix (Fin n) (Fin D) ℝ |
      SimpleSpec (Yᵀ * Y) (isHermitian_transpose_mul_self Y) rk}
      = ⋂ k : Fin (Fintype.card (Fin D)), ⋂ l : Fin (Fintype.card (Fin D)),
        {Y : Matrix (Fin n) (Fin D) ℝ |
          (k : ℕ) < rk → k ≠ l → gramEig Y (l : ℕ) ≠ gramEig Y (k : ℕ)} := by
    ext Y
    simp only [Set.mem_iInter, Set.mem_ofPred_eq]
    constructor
    · intro h k l hk hkl
      rw [gramEig_of_lt Y l.isLt, gramEig_of_lt Y k.isLt]
      have hl : (⟨(l : ℕ), l.isLt⟩ : Fin (Fintype.card (Fin D))) = l := Fin.eta _ _
      have hk' : (⟨(k : ℕ), k.isLt⟩ : Fin (Fintype.card (Fin D))) = k := Fin.eta _ _
      rw [hl, hk']
      exact h k l hk hkl
    · intro h k l hk hkl
      have hh := h k l hk hkl
      rw [gramEig_of_lt Y l.isLt, gramEig_of_lt Y k.isLt] at hh
      have hl : (⟨(l : ℕ), l.isLt⟩ : Fin (Fintype.card (Fin D))) = l := Fin.eta _ _
      have hk' : (⟨(k : ℕ), k.isLt⟩ : Fin (Fintype.card (Fin D))) = k := Fin.eta _ _
      rw [hl, hk'] at hh
      exact hh
  rw [hset]
  refine MeasurableSet.iInter fun k => MeasurableSet.iInter fun l => ?_
  by_cases hk : (k : ℕ) < rk
  · by_cases hkl : k ≠ l
    · have hset2 : {Y : Matrix (Fin n) (Fin D) ℝ |
          (k : ℕ) < rk → k ≠ l → gramEig Y (l : ℕ) ≠ gramEig Y (k : ℕ)}
          = (fun Y : Matrix (Fin n) (Fin D) ℝ =>
              gramEig Y (l : ℕ) - gramEig Y (k : ℕ)) ⁻¹' {(0 : ℝ)}ᶜ := by
        ext Y
        simp only [Set.mem_ofPred_eq, Set.mem_preimage, Set.mem_compl_iff,
          Set.mem_singleton_iff, sub_eq_zero]
        exact ⟨fun h => h hk hkl, fun h _ _ => h⟩
      rw [hset2]
      exact ((measurable_gramEig (l : ℕ)).sub (measurable_gramEig (k : ℕ)))
        (measurableSet_singleton (0 : ℝ)).compl
    · have : {Y : Matrix (Fin n) (Fin D) ℝ |
          (k : ℕ) < rk → k ≠ l → gramEig Y (l : ℕ) ≠ gramEig Y (k : ℕ)} = Set.univ := by
        ext Y
        simp only [Set.mem_ofPred_eq, Set.mem_univ, iff_true]
        intro _ hc
        exact absurd hc hkl
      rw [this]
      exact MeasurableSet.univ
  · have : {Y : Matrix (Fin n) (Fin D) ℝ |
        (k : ℕ) < rk → k ≠ l → gramEig Y (l : ℕ) ≠ gramEig Y (k : ℕ)} = Set.univ := by
      ext Y
      simp only [Set.mem_ofPred_eq, Set.mem_univ, iff_true]
      intro hc
      exact absurd hc hk
    rw [this]
    exact MeasurableSet.univ

/-- The same set for a matrix whose rows are indexed by a sum type, which is the shape
`YmatR` produces. The reindexing `Fin rr ⊕ Fin pN ≃ Fin (rr + pN)` leaves the Gram matrix
unchanged. -/
theorem measurableSet_simpleSpec_gram_sum {rr pN D : ℕ} (rk : ℕ) :
    MeasurableSet {Y : Matrix (Fin rr ⊕ Fin pN) (Fin D) ℝ |
      SimpleSpec (Yᵀ * Y) (isHermitian_transpose_mul_self' Y) rk} := by
  classical
  set e : Fin rr ⊕ Fin pN ≃ Fin (rr + pN) := finSumFinEquiv with hedef
  have hgram : ∀ Y : Matrix (Fin rr ⊕ Fin pN) (Fin D) ℝ,
      (Y.submatrix e.symm id)ᵀ * (Y.submatrix e.symm id) = Yᵀ * Y := by
    intro Y
    ext a b
    simp only [Matrix.mul_apply, Matrix.transpose_apply, Matrix.submatrix_apply, id_eq]
    exact Equiv.sum_comp e.symm fun k => Y k a * Y k b
  have hm : Measurable fun Y : Matrix (Fin rr ⊕ Fin pN) (Fin D) ℝ =>
      Y.submatrix e.symm id :=
    measurable_pi_lambda _ fun i => measurable_pi_lambda _ fun l =>
      (measurable_pi_apply l).comp (measurable_pi_apply (e.symm i))
  have hset : {Y : Matrix (Fin rr ⊕ Fin pN) (Fin D) ℝ |
      SimpleSpec (Yᵀ * Y) (isHermitian_transpose_mul_self' Y) rk}
      = (fun Y : Matrix (Fin rr ⊕ Fin pN) (Fin D) ℝ => Y.submatrix e.symm id) ⁻¹'
        {Z : Matrix (Fin (rr + pN)) (Fin D) ℝ |
          SimpleSpec (Zᵀ * Z) (isHermitian_transpose_mul_self Z) rk} := by
    ext Y
    simp only [Set.mem_ofPred_eq, Set.mem_preimage]
    exact ⟨fun h => simpleSpec_congr (hgram Y).symm _ _ h,
      fun h => simpleSpec_congr (hgram Y) _ _ h⟩
  rw [hset]
  exact hm (measurableSet_simpleSpec_gram rk)

/-- Joint measurability of `YmatR` in the two blocks. -/
theorem measurable_YmatR₂ {rr p D : ℕ} (t : ℝ) :
    Measurable fun ξ : Matrix (Fin rr) (Fin D) ℝ × Matrix (Fin p) (Fin D) ℝ =>
      YmatR ξ.1 t ξ.2 := by
  refine measurable_pi_lambda _ fun i => ?_
  cases i with
  | inl k =>
      have h : (fun ξ : Matrix (Fin rr) (Fin D) ℝ × Matrix (Fin p) (Fin D) ℝ =>
          YmatR ξ.1 t ξ.2 (Sum.inl k)) = fun ξ => ξ.1 k := by
        funext ξ; rfl
      rw [h]
      exact (measurable_pi_apply k).comp measurable_fst
  | inr i' =>
      refine measurable_pi_lambda _ fun l => ?_
      have h : (fun ξ : Matrix (Fin rr) (Fin D) ℝ × Matrix (Fin p) (Fin D) ℝ =>
          YmatR ξ.1 t ξ.2 (Sum.inr i') l) = fun ξ => t * ξ.2 i' l := by
        funext ξ; rfl
      rw [h]
      exact measurable_const.mul ((measurable_pi_apply l).comp
        ((measurable_pi_apply i').comp measurable_snd))

/-- Measurability of the `ℝ≥0∞`-valued integrand of the delocalization bound. The Gaussian
block is the second coordinate and the test direction reads the first. -/
theorem measurable_ofReal_overlapIdx₂ {rr pN D : ℕ} (a : ℕ)
    {g : Matrix (Fin rr) (Fin D) ℝ → EuclideanSpace ℝ (Fin D)} (hg : Measurable g) :
    Measurable fun ξ : Matrix (Fin rr) (Fin D) ℝ × Matrix (Fin pN) (Fin D) ℝ =>
      ENNReal.ofReal (overlapIdx ξ.2 a (g ξ.1)) := by
  have h1 : Measurable fun ξ : Matrix (Fin rr) (Fin D) ℝ × Matrix (Fin pN) (Fin D) ℝ =>
      (ξ.2, g ξ.1) := measurable_snd.prodMk (hg.comp measurable_fst)
  have h2 := (measurable_overlapIdx₂ (n := pN) (p := D) a).comp h1
  exact ENNReal.measurable_ofReal.comp h2

/-- `YmatR` with its rows reindexed by `Fin rr ⊕ Fin pN ≃ Fin (rr + pN)`. The reindexing
leaves the Gram matrix unchanged and puts the matrix in the shape `overlapIdx` takes. -/
noncomputable def ymatFlat {rr pN D : ℕ} (Y : Matrix (Fin rr) (Fin D) ℝ) (t : ℝ)
    (Bm : Matrix (Fin pN) (Fin D) ℝ) : Matrix (Fin (rr + pN)) (Fin D) ℝ :=
  (YmatR Y t Bm).submatrix finSumFinEquiv.symm id

theorem gram_ymatFlat {rr pN D : ℕ} (Y : Matrix (Fin rr) (Fin D) ℝ) (t : ℝ)
    (Bm : Matrix (Fin pN) (Fin D) ℝ) :
    (ymatFlat Y t Bm)ᵀ * ymatFlat Y t Bm = (YmatR Y t Bm)ᵀ * YmatR Y t Bm := by
  ext a b
  simp only [ymatFlat, Matrix.mul_apply, Matrix.transpose_apply, Matrix.submatrix_apply,
    id_eq]
  exact Equiv.sum_comp finSumFinEquiv.symm fun k => YmatR Y t Bm k a * YmatR Y t Bm k b

theorem measurable_ymatFlat₂ {rr pN D : ℕ} (t : ℝ) :
    Measurable fun ξ : Matrix (Fin rr) (Fin D) ℝ × Matrix (Fin pN) (Fin D) ℝ =>
      ymatFlat ξ.1 t ξ.2 :=
  (measurable_pi_lambda _ fun i => measurable_pi_lambda _ fun l =>
    (measurable_pi_apply l).comp
      (measurable_pi_apply (finSumFinEquiv.symm i))).comp (measurable_YmatR₂ t)

/-- Measurability of the `ℝ≥0∞`-valued integrand of the symmetry bound. The first coordinate
carries both the frame block of `YmatR` and the test direction. -/
theorem measurable_ofReal_overlapIdx_ymatFlat {rr pN D : ℕ} (k : ℕ) (t : ℝ)
    {Y : Matrix (Fin rr) (Fin D) ℝ → Matrix (Fin rr) (Fin D) ℝ} (hY : Measurable Y)
    {g : Matrix (Fin rr) (Fin D) ℝ → EuclideanSpace ℝ (Fin D)} (hg : Measurable g) :
    Measurable fun ξ : Matrix (Fin rr) (Fin D) ℝ × Matrix (Fin pN) (Fin D) ℝ =>
      ENNReal.ofReal (overlapIdx (ymatFlat (Y ξ.1) t ξ.2) k (g ξ.1)) := by
  have h1 : Measurable fun ξ : Matrix (Fin rr) (Fin D) ℝ × Matrix (Fin pN) (Fin D) ℝ =>
      (ymatFlat (Y ξ.1) t ξ.2, g ξ.1) :=
    ((measurable_ymatFlat₂ t).comp ((hY.comp measurable_fst).prodMk measurable_snd)).prodMk
      (hg.comp measurable_fst)
  have h2 := (measurable_overlapIdx₂ (n := rr + pN) (p := D) k).comp h1
  exact ENNReal.measurable_ofReal.comp h2

/-- A family of events of measure zero has measure tending to zero. -/
theorem tendsto_measure_zero_of_eventually_zero {Ω : ℕ → Type*}
    [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {s : ∀ N, Set (Ω N)}
    (h : ∀ᶠ N in atTop, μ N (s N) = 0) :
    Tendsto (fun N => μ N (s N)) atTop (𝓝 0) :=
  Tendsto.congr' (h.mono fun _ hN => hN.symm) tendsto_const_nhds

/-- Good sets whose complements sit eventually inside a vanishing family fill the space. -/
theorem tendsto_measure_one_of_bad_eventually {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} [∀ N, IsProbabilityMeasure (μ N)] {s t : ∀ N, Set (Ω N)}
    (hst : ∀ᶠ N in atTop, (t N)ᶜ ⊆ s N)
    (hs : Tendsto (fun N => μ N (s N)) atTop (𝓝 0)) :
    Tendsto (fun N => μ N (t N)) atTop (𝓝 1) := by
  have hone : Tendsto (fun _ : ℕ => (1 : ℝ≥0∞)) atTop (𝓝 1) := tendsto_const_nhds
  have hg : Tendsto (fun N => 1 - μ N (s N)) atTop (𝓝 1) := by
    have h1 : Tendsto (fun N => (1 : ℝ≥0∞) - μ N (s N)) atTop (𝓝 (1 - 0)) :=
      ENNReal.Tendsto.sub tendsto_const_nhds hs (Or.inl ENNReal.one_ne_top)
    simpa using h1
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' hg hone ?_
    (Filter.Eventually.of_forall fun _ => prob_le_one)
  filter_upwards [hst] with N hN
  refine tsub_le_iff_right.mpr ?_
  calc (1 : ℝ≥0∞) = μ N Set.univ := measure_univ.symm
    _ = μ N (t N ∪ (t N)ᶜ) := by rw [Set.union_compl_self]
    _ ≤ μ N (t N) + μ N ((t N)ᶜ) := measure_union_le _ _
    _ ≤ μ N (t N) + μ N (s N) := add_le_add le_rfl (measure_mono hN)

end Helpers

end EdgeGlueR

/-! ### 2. The split with the pair law kept -/

namespace RankRStack

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-- `resolventLimitsR_of_gaussian` (`RankR/RMT/Forms.lean`) with the Gaussian block `B` and its
joint law with the frame coefficients kept in the conclusion. The body is the body of that
theorem; only the last conjunct is new. Task U7c needs `B` because the delocalization bound and
the symmetry bound are conditional statements on the pair `(Uᵀ Zu, B)`. -/
theorem resolventLimitsR_of_gaussian_pair [∀ N, IsProbabilityMeasure (μ N)]
    (s : RankRStack μ ns d r) (hG : s.GaussianNoise) {c : ℝ} (hc : 0 < c)
    (hdtop : Tendsto d atTop atTop)
    (hns : Tendsto (fun N => (ns N : ℝ) / d N) atTop (𝓝 c))
    {pp : ℕ → ℕ} (hn : ∀ N, ns N = r + pp N) (hp : ∀ N, 0 < pp N) :
    ∃ U : (N : ℕ) → Matrix (Fin (ns N)) (Fin r) ℝ,
      (∀ N, (U N)ᵀ * U N = 1) ∧
      (∀ N ω, s.gram N ω
        = s.rankRW0 N ω (U N) + s.qmatR N ω (U N) * (s.qmatR N ω (U N))ᵀ) ∧
      ResolventLimitsR μ (fun N ω => s.rankRW0 N ω (U N))
        (fun N ω => s.isHermitian_rankRW0 N ω (U N)) s.spikeVec
        (fun k N ω => fun l => ((s.E N ω)ᵀ * U N) l k) c ∧
      (∀ N, ∃ B : Ω N → Matrix (Fin (pp N)) (Fin (d N)) ℝ,
        (∀ ω, s.rankRW0 N ω (U N) = ((d N : ℝ))⁻¹ • ((B ω)ᵀ * B ω)) ∧
        HasLaw (fun ω => ((U N)ᵀ * s.Zu N ω, B ω))
          ((gaussianMatrix r (d N)).prod (gaussianMatrix (pp N) (d N))) (μ N)) := by
  classical
  choose U B hU hsig hgram hlaw using fun N => s.exists_rankR_split hG N (hn N)
  have hd0 : ∀ N, 0 < d N := s.hd
  have hr : ∀ N, r ≤ ns N := fun N => by
    rw [hn N]
    exact Nat.le_add_right r (pp N)
  have hppN : ∀ N, ns N - r = pp N := fun N => by
    rw [hn N]
    omega
  have hcN : Tendsto (fun N => (pp N : ℝ) / d N) atTop (𝓝 c) := by
    have h := tendsto_sub_ratio (nf := ns) r hr hdtop hns
    simpa only [hppN] using h
  set ZZ : ∀ N, Ω N → Matrix (Fin r) (Fin (d N)) ℝ × Matrix (Fin (pp N)) (Fin (d N)) ℝ :=
    fun N ω => ((U N)ᵀ * s.Zu N ω, B N ω) with hZZdef
  have hZZlaw : ∀ N, HasLaw (ZZ N)
      ((gaussianMatrix r (d N)).prod (gaussianMatrix (pp N) (d N))) (μ N) := hlaw
  have hBsnd : ∀ N, HasLaw (B N) (gaussianMatrix (pp N) (d N)) (μ N) :=
    fun N => measurePreserving_snd.fun_comp_hasLaw (hlaw N)
  have hW0eq : ∀ N ω, R2.W0 (ZZ N ω).2 = s.rankRW0 N ω (U N) := fun N ω =>
    (s.rankRW0_eq_smul_block N ω (hU N) (hsig N) (hgram N ω)).symm
  have hgeq : ∀ (k : Fin r) (N : ℕ) (ω : Ω N),
      R2.gOf ((ZZ N ω).1 k) = fun l => ((s.E N ω)ᵀ * U N) l k :=
    fun k N ω => s.gOf_row_eq_col N ω (U N) k
  have hedge1 : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
        ≤ bulkEdge c + ε}) atTop (𝓝 1) := by
    intro ε hε
    exact R3.tendsto_measure_lamMax_le μ hc B (fun N ω => s.rankRW0 N ω (U N))
      (fun N ω => s.rankRW0_eq_smul_block N ω (hU N) (hsig N) (hgram N ω))
      (fun N ω => s.isHermitian_rankRW0 N ω (U N)) hBsnd hp hd0 hdtop hcN hε
  have hedgeC : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
        ≤ bulkEdge c + ε}ᶜ) atTop (𝓝 0) := fun ε hε =>
    tendsto_measure_compl_zero
      (fun N => (s.measurableSet_lamMax_rankRW0_le N (U N) _).nullMeasurableSet)
      (hedge1 ε hε)
  refine ⟨U, hU, fun N ω => s.gram_eq_rankRW0_add N ω (hU N) (hsig N), ?_, ?_⟩
  · exact resolventLimitsR_of_pairLaw hc hdtop hp hcN (fun N ω => s.rankRW0 N ω (U N))
      (fun N ω => s.isHermitian_rankRW0 N ω (U N)) s.spikeVec
      (fun k N ω => fun l => ((s.E N ω)ᵀ * U N) l k) ZZ hZZlaw hW0eq hgeq
      (fun a b N => s.inner_spikeVec N a b) hedge1 hedgeC
  · intro N
    exact ⟨B N, fun ω => s.rankRW0_eq_smul_block N ω (hU N) (hsig N) (hgram N ω), hlaw N⟩

/-! `measurable_X` and `measurableSet_eigenvalues₀_le` used to stand here, in a second copy
next to the one of `RankR/RMT/AlignOutG.lean`. The second cleanup pass (2026-09-02) moved one
copy of each to `RankR/RMT/EdgeR.lean`, which both files import. -/

/-- The bad-set form of `tendsto_measure_eigenvalues₀_le` (`RankR/RMT/EdgeR.lean`). -/
theorem tendsto_measure_eigenvalues₀_gt [∀ N, IsProbabilityMeasure (μ N)]
    (s : RankRStack μ ns d r) {U : (N : ℕ) → Matrix (Fin (ns N)) (Fin r) ℝ} {c : ℝ}
    (hc : 0 < c)
    (h : ResolventLimitsR μ (fun N ω => s.rankRW0 N ω (U N))
      (fun N ω => s.isHermitian_rankRW0 N ω (U N)) s.spikeVec
      (fun k N ω => fun l => ((s.E N ω)ᵀ * U N) l k) c)
    (hgram : ∀ N ω, s.gram N ω
      = s.rankRW0 N ω (U N) + s.qmatR N ω (U N) * (s.qmatR N ω (U N))ᵀ)
    {t u : ℕ} (e : Fin t ⊕ Fin u ≃ Fin r)
    (hsub : ∀ a : Fin t, s.coreEig (e (Sum.inl a)) ^ 2 ≤ c) :
    ∀ ε > 0, Tendsto (fun N => μ N
        {ω | ∀ k : Fin (Fintype.card (Fin (d N))), u ≤ (k : ℕ) →
          (s.isHermitian_gram N ω).eigenvalues₀ k ≤ bulkEdge c + ε}ᶜ) atTop (𝓝 0) := by
  intro ε hε
  exact tendsto_measure_compl_zero
    (fun N => (s.measurableSet_eigenvalues₀_le N u (bulkEdge c + ε)).nullMeasurableSet)
    (s.tendsto_measure_eigenvalues₀_le hc h hgram e hsub ε hε)


/-! ### 3. The objects of the conditioning, as functions of the frame coefficients -/

section Conditioning

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-- `Q` written as a function of the frame coefficients `α = Uᵀ Zu` alone. -/
noncomputable def qmatOf (s : RankRStack μ ns d r) (N : ℕ)
    (α : Matrix (Fin r) (Fin (d N)) ℝ) : Matrix (Fin (d N)) (Fin r) ℝ :=
  s.spikeMat N * Matrix.diagonal (fun j => Real.sqrt (s.coreEig j))
    + (Real.sqrt (d N))⁻¹ • αᵀ

theorem qmatOf_eq (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (ns N)) (Fin r) ℝ) : s.qmatOf N (Uᵀ * s.Zu N ω) = s.qmatR N ω U := by
  rw [qmatOf, qmatR, s.E_transpose_mul N ω U]

/-- Column `l` of `Q`, as a Euclidean vector. -/
noncomputable def colOf (s : RankRStack μ ns d r) (N : ℕ)
    (α : Matrix (Fin r) (Fin (d N)) ℝ) (l : Fin r) : EuclideanSpace ℝ (Fin (d N)) :=
  WithLp.toLp 2 fun i => s.qmatOf N α i l

/-- The column truncated at level `C`. See the header note on the truncation. -/
noncomputable def colTruncOf (s : RankRStack μ ns d r) (N : ℕ)
    (α : Matrix (Fin r) (Fin (d N)) ℝ) (l : Fin r) (C : ℝ) : EuclideanSpace ℝ (Fin (d N)) :=
  if ‖s.colOf N α l‖ ^ 2 ≤ C then s.colOf N α l else 0

theorem norm_colTruncOf_sq_le (s : RankRStack μ ns d r) (N : ℕ)
    (α : Matrix (Fin r) (Fin (d N)) ℝ) (l : Fin r) {C : ℝ} (hC : 0 ≤ C) :
    ‖s.colTruncOf N α l C‖ ^ 2 ≤ C := by
  rw [colTruncOf]
  split
  · assumption
  · simpa using hC

theorem colTruncOf_eq_colOf (s : RankRStack μ ns d r) (N : ℕ)
    (α : Matrix (Fin r) (Fin (d N)) ℝ) (l : Fin r) {C : ℝ} (h : ‖s.colOf N α l‖ ^ 2 ≤ C) :
    s.colTruncOf N α l C = s.colOf N α l := by
  rw [colTruncOf, if_pos h]

/-- The residual direction `v - Q a` of the decomposition, as a function of `α`, with the junk
value `0` where `QᵀQ` is singular. -/
noncomputable def rvecOf (s : RankRStack μ ns d r) (N : ℕ) (j : Fin r)
    (α : Matrix (Fin r) (Fin (d N)) ℝ) : EuclideanSpace ℝ (Fin (d N)) :=
  Set.indicator {β : Matrix (Fin r) (Fin (d N)) ℝ |
      IsUnit ((s.qmatOf N β)ᵀ * s.qmatOf N β).det}
    (fun β => WithLp.toLp 2
      (DecompR.rvecOfR (s.qmatOf N β) (WithLp.ofLp (s.spikeVec j N)))) α

theorem rvecOf_of_isUnit (s : RankRStack μ ns d r) (N : ℕ) (j : Fin r)
    {α : Matrix (Fin r) (Fin (d N)) ℝ} (h : IsUnit ((s.qmatOf N α)ᵀ * s.qmatOf N α).det) :
    s.rvecOf N j α
      = WithLp.toLp 2 (DecompR.rvecOfR (s.qmatOf N α) (WithLp.ofLp (s.spikeVec j N))) := by
  rw [rvecOf]
  exact Set.indicator_of_mem (a := α) h _

theorem rvecOf_of_not_isUnit (s : RankRStack μ ns d r) (N : ℕ) (j : Fin r)
    {α : Matrix (Fin r) (Fin (d N)) ℝ} (h : ¬ IsUnit ((s.qmatOf N α)ᵀ * s.qmatOf N α).det) :
    s.rvecOf N j α = 0 := by
  rw [rvecOf]
  exact Set.indicator_of_notMem (a := α) h _

/-! #### Measurability -/

theorem measurable_qmatOf_apply (s : RankRStack μ ns d r) (N : ℕ) (i : Fin (d N)) (l : Fin r) :
    Measurable fun α : Matrix (Fin r) (Fin (d N)) ℝ => s.qmatOf N α i l := by
  have h : (fun α : Matrix (Fin r) (Fin (d N)) ℝ => s.qmatOf N α i l)
      = fun α => (s.spikeMat N * Matrix.diagonal fun j => Real.sqrt (s.coreEig j)) i l
        + (Real.sqrt (d N))⁻¹ * α l i := by
    funext α
    rfl
  rw [h]
  exact measurable_const.add (measurable_const.mul
    ((measurable_pi_apply i).comp (measurable_pi_apply l)))

theorem measurable_qmatOf (s : RankRStack μ ns d r) (N : ℕ) :
    Measurable fun α : Matrix (Fin r) (Fin (d N)) ℝ => s.qmatOf N α :=
  measurable_pi_lambda _ fun i => measurable_pi_lambda _ fun l =>
    s.measurable_qmatOf_apply N i l

theorem measurable_colOf (s : RankRStack μ ns d r) (N : ℕ) (l : Fin r) :
    Measurable fun α : Matrix (Fin r) (Fin (d N)) ℝ => s.colOf N α l :=
  (WithLp.measurable_toLp 2 (Fin (d N) → ℝ)).comp
    (measurable_pi_lambda _ fun i => s.measurable_qmatOf_apply N i l)

theorem measurable_colTruncOf (s : RankRStack μ ns d r) (N : ℕ) (l : Fin r) (C : ℝ) :
    Measurable fun α : Matrix (Fin r) (Fin (d N)) ℝ => s.colTruncOf N α l C := by
  have hnorm : Measurable fun α : Matrix (Fin r) (Fin (d N)) ℝ => ‖s.colOf N α l‖ ^ 2 :=
    (measurable_norm.comp (s.measurable_colOf N l)).pow_const 2
  have hset : MeasurableSet {α : Matrix (Fin r) (Fin (d N)) ℝ | ‖s.colOf N α l‖ ^ 2 ≤ C} :=
    hnorm measurableSet_Iic
  exact Measurable.ite hset (s.measurable_colOf N l) measurable_const

theorem measurable_det_qmatOf (s : RankRStack μ ns d r) (N : ℕ) :
    Measurable fun α : Matrix (Fin r) (Fin (d N)) ℝ =>
      ((s.qmatOf N α)ᵀ * s.qmatOf N α).det := by
  refine EdgeGlueR.measurable_det_of_entries fun i j => ?_
  have h : (fun α : Matrix (Fin r) (Fin (d N)) ℝ =>
      ((s.qmatOf N α)ᵀ * s.qmatOf N α) i j)
      = fun α => ∑ k, s.qmatOf N α k i * s.qmatOf N α k j := by
    funext α
    simp only [Matrix.mul_apply, Matrix.transpose_apply]
  rw [h]
  exact Finset.measurable_sum _ fun k _ =>
    (s.measurable_qmatOf_apply N k i).mul (s.measurable_qmatOf_apply N k j)

theorem measurableSet_isUnit_det_qmatOf (s : RankRStack μ ns d r) (N : ℕ) :
    MeasurableSet {α : Matrix (Fin r) (Fin (d N)) ℝ |
      IsUnit ((s.qmatOf N α)ᵀ * s.qmatOf N α).det} := by
  have hset : {α : Matrix (Fin r) (Fin (d N)) ℝ |
      IsUnit ((s.qmatOf N α)ᵀ * s.qmatOf N α).det}
      = (fun α => ((s.qmatOf N α)ᵀ * s.qmatOf N α).det) ⁻¹' {(0 : ℝ)}ᶜ := by
    ext α
    simp only [Set.mem_ofPred_eq, Set.mem_preimage, Set.mem_compl_iff, Set.mem_singleton_iff]
    exact isUnit_iff_ne_zero
  rw [hset]
  exact (s.measurable_det_qmatOf N) (measurableSet_singleton (0 : ℝ)).compl

theorem measurable_rvecOf (s : RankRStack μ ns d r) (N : ℕ) (j : Fin r) :
    Measurable fun α : Matrix (Fin r) (Fin (d N)) ℝ => s.rvecOf N j α := by
  have hf : Measurable fun α : Matrix (Fin r) (Fin (d N)) ℝ => (WithLp.toLp 2
      (DecompR.rvecOfR (s.qmatOf N α) (WithLp.ofLp (s.spikeVec j N)))
        : EuclideanSpace ℝ (Fin (d N))) :=
    (WithLp.measurable_toLp 2 (Fin (d N) → ℝ)).comp
      (DecompR.measurable_rvecOfR.comp ((s.measurable_qmatOf N).prodMk measurable_const))
  exact hf.indicator (s.measurableSet_isUnit_det_qmatOf N)


end Conditioning

/-! ### 4. The two almost sure simplicity events (G2 and G3) -/

section Simplicity

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-- The scaling of the noise, `t = d^{-1/2}`, is nonzero. -/
theorem sqrt_inv_ne_zero (s : RankRStack μ ns d r) (N : ℕ) :
    (Real.sqrt (d N))⁻¹ ≠ 0 := by
  have hdpos : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast s.hd N
  have h : Real.sqrt (d N) ≠ 0 := ne_of_gt (Real.sqrt_pos.mpr hdpos)
  exact inv_ne_zero h

/-- Measurability of the affine simplicity predicate. -/
theorem measurable_simpleSpec_affine {q D : ℕ} (A : Matrix (Fin q) (Fin D) ℝ) (t : ℝ)
    (rk : ℕ) :
    Measurable fun Z : Matrix (Fin q) (Fin D) ℝ =>
      SimpleSpec ((A + t • Z)ᵀ * (A + t • Z))
        (isHermitian_transpose_mul_self (A + t • Z)) rk := by
  rw [← measurableSet_setOfPred]
  exact measurable_affine_matrix A t (EdgeGlueR.measurableSet_simpleSpec_gram rk)

/-- **(G2)** The block `W₀` of the split has an almost surely simple top-`r` spectrum. -/
theorem ae_simpleSpec_rankRW0 (s : RankRStack μ ns d r) (N : ℕ)
    {U : Matrix (Fin (ns N)) (Fin r) ℝ} {pN : ℕ}
    {B : Ω N → Matrix (Fin pN) (Fin (d N)) ℝ}
    (hW0 : ∀ ω, s.rankRW0 N ω U = ((d N : ℝ))⁻¹ • ((B ω)ᵀ * B ω))
    (hBlaw : HasLaw B (gaussianMatrix pN (d N)) (μ N))
    (hpN : 0 < pN) (hrk : r ≤ min pN (d N)) :
    ∀ᵐ ω ∂(μ N), SimpleSpec (s.rankRW0 N ω U) (s.isHermitian_rankRW0 N ω U) r := by
  classical
  have ht : (Real.sqrt (d N))⁻¹ ≠ 0 := s.sqrt_inv_ne_zero N
  have hsq : (Real.sqrt (d N))⁻¹ * (Real.sqrt (d N))⁻¹ = ((d N : ℝ))⁻¹ :=
    R2.sqrt_inv_mul_sqrt_inv (d N)
  have hae := simpleSpec_ae_affine hpN (s.hd N) (0 : Matrix (Fin pN) (Fin (d N)) ℝ) ht r hrk
  have htrans := (hBlaw.ae_iff
    (measurable_simpleSpec_affine (0 : Matrix (Fin pN) (Fin (d N)) ℝ)
      (Real.sqrt (d N))⁻¹ r)).mpr hae
  filter_upwards [htrans] with ω hω
  refine EdgeGlueR.simpleSpec_congr ?_ _ _ hω
  rw [zero_add, Matrix.transpose_smul, Matrix.smul_mul, Matrix.mul_smul, smul_smul, hsq]
  exact (hW0 ω).symm

/-- **(G3)** The Gram matrix of the table has an almost surely simple top-`r` spectrum. -/
theorem ae_simpleSpec_gram (s : RankRStack μ ns d r) (hG : s.GaussianNoise) (N : ℕ)
    (hns0 : 0 < ns N) (hrk : r ≤ min (ns N) (d N)) :
    ∀ᵐ ω ∂(μ N), SimpleSpec (s.gram N ω) (s.isHermitian_gram N ω) r := by
  classical
  have ht : (Real.sqrt (d N))⁻¹ ≠ 0 := s.sqrt_inv_ne_zero N
  have hae := simpleSpec_ae_affine hns0 (s.hd N) (s.signalPart N) ht r hrk
  have htrans := ((hG N).ae_iff
    (measurable_simpleSpec_affine (s.signalPart N) (Real.sqrt (d N))⁻¹ r)).mpr hae
  filter_upwards [htrans] with ω hω
  refine EdgeGlueR.simpleSpec_congr ?_ _ _ hω
  have hX : s.signalPart N + (Real.sqrt (d N))⁻¹ • s.Zu N ω = s.X N ω := by
    rw [s.X_eq N ω, s.hE N ω]
  rw [hX]
  rfl

/-- The Gram matrix of the table, written through the `YmatR` of item Sym. -/
theorem gram_eq_YmatR (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    {U : Matrix (Fin (ns N)) (Fin r) ℝ} {pN : ℕ}
    {B : Ω N → Matrix (Fin pN) (Fin (d N)) ℝ}
    (hgram : s.gram N ω = s.rankRW0 N ω U + s.qmatR N ω U * (s.qmatR N ω U)ᵀ)
    (hW0 : s.rankRW0 N ω U = ((d N : ℝ))⁻¹ • ((B ω)ᵀ * B ω)) :
    s.gram N ω
      = (YmatR (s.qmatOf N (Uᵀ * s.Zu N ω))ᵀ (Real.sqrt (d N))⁻¹ (B ω))ᵀ
        * YmatR (s.qmatOf N (Uᵀ * s.Zu N ω))ᵀ (Real.sqrt (d N))⁻¹ (B ω) := by
  have hsq : ((Real.sqrt (d N))⁻¹ : ℝ) ^ 2 = ((d N : ℝ))⁻¹ := by
    rw [sq]
    exact R2.sqrt_inv_mul_sqrt_inv (d N)
  rw [gram_YmatR, hsq, Matrix.transpose_transpose, s.qmatOf_eq N ω U, hgram, hW0, add_comm]

/-- Joint measurability of the `YmatR` of the split in the two Gaussian blocks. -/
theorem measurable_YmatR_qmatOf (s : RankRStack μ ns d r) (N : ℕ) {pN : ℕ} :
    Measurable fun ξ : Matrix (Fin r) (Fin (d N)) ℝ × Matrix (Fin pN) (Fin (d N)) ℝ =>
      YmatR (s.qmatOf N ξ.1)ᵀ (Real.sqrt (d N))⁻¹ ξ.2 := by
  have hQ : Measurable fun ξ : Matrix (Fin r) (Fin (d N)) ℝ ×
      Matrix (Fin pN) (Fin (d N)) ℝ => (s.qmatOf N ξ.1)ᵀ :=
    (measurable_pi_lambda _ fun l => measurable_pi_lambda _ fun i =>
      (s.measurable_qmatOf_apply N i l).comp measurable_fst)
  exact (EdgeGlueR.measurable_YmatR₂ (Real.sqrt (d N))⁻¹).comp (hQ.prodMk measurable_snd)

/-- **(G3) on the product space.** The fiberwise simplicity hypothesis of
`lintegral_normSq_specProjTop_YmatR_le`. -/
theorem ae_ae_simpleSpec_YmatR (s : RankRStack μ ns d r) (hG : s.GaussianNoise) (N : ℕ)
    {U : Matrix (Fin (ns N)) (Fin r) ℝ} {pN : ℕ}
    {B : Ω N → Matrix (Fin pN) (Fin (d N)) ℝ}
    (hgram : ∀ ω, s.gram N ω = s.rankRW0 N ω U + s.qmatR N ω U * (s.qmatR N ω U)ᵀ)
    (hW0 : ∀ ω, s.rankRW0 N ω U = ((d N : ℝ))⁻¹ • ((B ω)ᵀ * B ω))
    (hlaw : HasLaw (fun ω => (Uᵀ * s.Zu N ω, B ω))
      ((gaussianMatrix r (d N)).prod (gaussianMatrix pN (d N))) (μ N))
    (hns0 : 0 < ns N) (hrk : r ≤ min (ns N) (d N)) :
    ∀ᵐ α ∂(gaussianMatrix r (d N)), ∀ᵐ Bm ∂(gaussianMatrix pN (d N)),
      SimpleSpec ((YmatR (s.qmatOf N α)ᵀ (Real.sqrt (d N))⁻¹ Bm)ᵀ
          * YmatR (s.qmatOf N α)ᵀ (Real.sqrt (d N))⁻¹ Bm)
        (isHermitian_transpose_mul_self' _) r := by
  classical
  have hpred : Measurable fun ξ : Matrix (Fin r) (Fin (d N)) ℝ ×
      Matrix (Fin pN) (Fin (d N)) ℝ =>
      SimpleSpec ((YmatR (s.qmatOf N ξ.1)ᵀ (Real.sqrt (d N))⁻¹ ξ.2)ᵀ
          * YmatR (s.qmatOf N ξ.1)ᵀ (Real.sqrt (d N))⁻¹ ξ.2)
        (isHermitian_transpose_mul_self' _) r := by
    rw [← measurableSet_setOfPred]
    exact (s.measurable_YmatR_qmatOf N) (EdgeGlueR.measurableSet_simpleSpec_gram_sum r)
  have hmu : ∀ᵐ ω ∂(μ N),
      SimpleSpec ((YmatR (s.qmatOf N (Uᵀ * s.Zu N ω))ᵀ (Real.sqrt (d N))⁻¹ (B ω))ᵀ
          * YmatR (s.qmatOf N (Uᵀ * s.Zu N ω))ᵀ (Real.sqrt (d N))⁻¹ (B ω))
        (isHermitian_transpose_mul_self' _) r := by
    filter_upwards [s.ae_simpleSpec_gram hG N hns0 hrk] with ω hω
    exact EdgeGlueR.simpleSpec_congr (s.gram_eq_YmatR N ω (hgram ω) (hW0 ω)) _ _ hω
  exact Measure.ae_ae_of_ae_prod ((hlaw.ae_iff hpred).mp hmu)


end Simplicity

/-! ### 5. The delocalization event (G4) -/

section Deloc

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-- The projector at a sorted index of `W₀` is the projector of the unscaled Wishart block.
The scaling bridge `EdgeGlueDetR.specProjIdx_smul` removes the factor `d⁻¹`. -/
theorem specProjIdx_rankRW0_eq (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (ns N)) (Fin r) ℝ) {pN : ℕ} {Bm : Matrix (Fin pN) (Fin (d N)) ℝ}
    (hW0 : s.rankRW0 N ω U = ((d N : ℝ))⁻¹ • (Bmᵀ * Bm)) (a : ℕ) :
    specProjIdx (s.rankRW0 N ω U) (s.isHermitian_rankRW0 N ω U) a
      = specProjIdx (Bmᵀ * Bm) (isHermitian_transpose_mul_self Bm) a := by
  have h' : (((d N : ℝ))⁻¹ • (Bmᵀ * Bm)).IsHermitian := by
    rw [← hW0]
    exact s.isHermitian_rankRW0 N ω U
  have hdpos : (0 : ℝ) < ((d N : ℝ))⁻¹ := by
    have : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast s.hd N
    exact inv_pos.mpr this
  rw [specProjIdx_congr hW0 (s.isHermitian_rankRW0 N ω U) h' a]
  exact EdgeGlueDetR.specProjIdx_smul (isHermitian_transpose_mul_self Bm) h' hdpos a

/-- **(G4), Markov form.** The double sum of the squared projections of the truncated columns
of `Q` on the top-`r` eigenvectors of `W₀` exceeds `κ` with probability at most
`r (∑ C) / (d κ)`. -/
theorem measure_kappa_ge_le (s : RankRStack μ ns d r) (N : ℕ)
    {U : Matrix (Fin (ns N)) (Fin r) ℝ} {pN : ℕ}
    {B : Ω N → Matrix (Fin pN) (Fin (d N)) ℝ}
    (hW0 : ∀ ω, s.rankRW0 N ω U = ((d N : ℝ))⁻¹ • ((B ω)ᵀ * B ω))
    (hlaw : HasLaw (fun ω => (Uᵀ * s.Zu N ω, B ω))
      ((gaussianMatrix r (d N)).prod (gaussianMatrix pN (d N))) (μ N))
    (hrp : r ≤ min pN (d N)) {CC : Fin r → ℝ} (hCC : ∀ l, 0 ≤ CC l) {κ : ℝ} (hκ : 0 < κ) :
    μ N {ω | κ ≤ ∑ a ∈ Finset.range r, ∑ l,
        ‖specProjIdx (s.rankRW0 N ω U) (s.isHermitian_rankRW0 N ω U) a
          (s.colTruncOf N (Uᵀ * s.Zu N ω) l (CC l))‖ ^ 2}
      ≤ ENNReal.ofReal ((r : ℝ) * (∑ l, CC l) / ((d N : ℝ) * κ)) := by
  classical
  have hdpos : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast s.hd N
  set F : Matrix (Fin r) (Fin (d N)) ℝ × Matrix (Fin pN) (Fin (d N)) ℝ → ℝ≥0∞ :=
    fun ξ => ∑ a ∈ Finset.range r, ∑ l,
      ENNReal.ofReal (overlapIdx ξ.2 a (s.colTruncOf N ξ.1 l (CC l))) with hFdef
  -- measurability of `F`
  have hterm : ∀ (a : ℕ) (l : Fin r), Measurable fun ξ : Matrix (Fin r) (Fin (d N)) ℝ ×
      Matrix (Fin pN) (Fin (d N)) ℝ => ENNReal.ofReal
        (overlapIdx ξ.2 a (s.colTruncOf N ξ.1 l (CC l))) := by
    intro a l
    exact EdgeGlueR.measurable_ofReal_overlapIdx₂ (g := fun α => s.colTruncOf N α l (CC l)) a
      (s.measurable_colTruncOf N l (CC l))
  have hFm : Measurable F := by
    rw [hFdef]
    exact Finset.measurable_sum _ fun a _ => Finset.measurable_sum _ fun l _ => hterm a l
  -- the pointwise identification
  have hpt : ∀ ω, ENNReal.ofReal (∑ a ∈ Finset.range r, ∑ l,
      ‖specProjIdx (s.rankRW0 N ω U) (s.isHermitian_rankRW0 N ω U) a
        (s.colTruncOf N (Uᵀ * s.Zu N ω) l (CC l))‖ ^ 2)
      = F (Uᵀ * s.Zu N ω, B ω) := by
    intro ω
    rw [hFdef]
    rw [ENNReal.ofReal_sum_of_nonneg
      (fun a _ => Finset.sum_nonneg fun l _ => sq_nonneg _)]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [ENNReal.ofReal_sum_of_nonneg (fun l _ => sq_nonneg _)]
    refine Finset.sum_congr rfl fun l _ => ?_
    congr 1
    rw [overlapIdx, ← s.specProjIdx_rankRW0_eq N ω U (hW0 ω) a]
  -- Markov
  have hmeas : AEMeasurable (fun ω => ENNReal.ofReal (∑ a ∈ Finset.range r, ∑ l,
      ‖specProjIdx (s.rankRW0 N ω U) (s.isHermitian_rankRW0 N ω U) a
        (s.colTruncOf N (Uᵀ * s.Zu N ω) l (CC l))‖ ^ 2)) (μ N) := by
    refine AEMeasurable.congr ?_ (Filter.Eventually.of_forall fun ω => (hpt ω).symm)
    exact hFm.aemeasurable.comp_aemeasurable hlaw.aemeasurable
  have hset : {ω | κ ≤ ∑ a ∈ Finset.range r, ∑ l,
        ‖specProjIdx (s.rankRW0 N ω U) (s.isHermitian_rankRW0 N ω U) a
          (s.colTruncOf N (Uᵀ * s.Zu N ω) l (CC l))‖ ^ 2}
      = {ω | ENNReal.ofReal κ ≤ ENNReal.ofReal (∑ a ∈ Finset.range r, ∑ l,
          ‖specProjIdx (s.rankRW0 N ω U) (s.isHermitian_rankRW0 N ω U) a
            (s.colTruncOf N (Uᵀ * s.Zu N ω) l (CC l))‖ ^ 2)} := by
    ext ω
    simp only [Set.mem_ofPred_eq]
    rw [ENNReal.ofReal_le_ofReal_iff
      (Finset.sum_nonneg fun a _ => Finset.sum_nonneg fun l _ => sq_nonneg _)]
  have hκne : ENNReal.ofReal κ ≠ 0 := by simp [ENNReal.ofReal_eq_zero, not_le, hκ]
  rw [hset]
  refine (meas_ge_le_lintegral_div hmeas hκne ENNReal.ofReal_ne_top).trans ?_
  -- the integral over the product measure
  have hint : ∫⁻ ω, ENNReal.ofReal (∑ a ∈ Finset.range r, ∑ l,
      ‖specProjIdx (s.rankRW0 N ω U) (s.isHermitian_rankRW0 N ω U) a
        (s.colTruncOf N (Uᵀ * s.Zu N ω) l (CC l))‖ ^ 2) ∂(μ N)
      = ∫⁻ ξ, F ξ ∂((gaussianMatrix r (d N)).prod (gaussianMatrix pN (d N))) := by
    rw [lintegral_congr hpt]
    exact hlaw.lintegral_comp hFm.aemeasurable
  have hbound : ∫⁻ ξ, F ξ ∂((gaussianMatrix r (d N)).prod (gaussianMatrix pN (d N)))
      ≤ ENNReal.ofReal ((r : ℝ) * (∑ l, CC l) / (d N : ℝ)) := by
    rw [hFdef]
    rw [lintegral_finsetSum _ (fun a _ => Finset.measurable_sum _ fun l _ => hterm a l)]
    have hin : ∀ a ∈ Finset.range r,
        ∫⁻ ξ, ∑ l, ENNReal.ofReal (overlapIdx ξ.2 a (s.colTruncOf N ξ.1 l (CC l)))
            ∂((gaussianMatrix r (d N)).prod (gaussianMatrix pN (d N)))
          ≤ ENNReal.ofReal ((∑ l, CC l) / (d N : ℝ)) := by
      intro a ha
      rw [lintegral_finsetSum _ (fun l _ => hterm a l)]
      have hone : ∀ l : Fin r,
          ∫⁻ ξ, ENNReal.ofReal (overlapIdx ξ.2 a (s.colTruncOf N ξ.1 l (CC l)))
              ∂((gaussianMatrix r (d N)).prod (gaussianMatrix pN (d N)))
            ≤ ENNReal.ofReal (CC l / (d N : ℝ)) := by
        intro l
        have hak : a < min pN (d N) := lt_of_lt_of_le (Finset.mem_range.mp ha) hrp
        set q : Matrix (Fin r) (Fin (d N)) ℝ → EuclideanSpace ℝ (Fin (d N)) :=
          fun α => s.colTruncOf N α l (CC l) with hqdef
        set G : Matrix (Fin pN) (Fin (d N)) ℝ × Matrix (Fin r) (Fin (d N)) ℝ → ℝ≥0∞ :=
          fun ζ => ENNReal.ofReal (‖specProjIdx (ζ.1ᵀ * ζ.1)
            (isHermitian_transpose_mul_self ζ.1) a (q ζ.2)‖ ^ 2) with hGdef
        have hswap := lintegral_prod_swap (μ := gaussianMatrix pN (d N))
          (ν := gaussianMatrix r (d N)) G
        have hrw : ∫⁻ ξ, ENNReal.ofReal (overlapIdx ξ.2 a (s.colTruncOf N ξ.1 l (CC l)))
            ∂((gaussianMatrix r (d N)).prod (gaussianMatrix pN (d N)))
            = ∫⁻ ζ, G ζ ∂((gaussianMatrix pN (d N)).prod (gaussianMatrix r (d N))) := by
          rw [← hswap]
          rfl
        rw [hrw]
        refine (lintegral_prod_normSq_specProjIdx_le hak (gaussianMatrix r (d N)) q).trans ?_
        have hqb : ∫⁻ α, ENNReal.ofReal (‖q α‖ ^ 2) ∂(gaussianMatrix r (d N))
            ≤ ENNReal.ofReal (CC l) := by
          calc ∫⁻ α, ENNReal.ofReal (‖q α‖ ^ 2) ∂(gaussianMatrix r (d N))
              ≤ ∫⁻ _, ENNReal.ofReal (CC l) ∂(gaussianMatrix r (d N)) :=
                lintegral_mono fun α => ENNReal.ofReal_le_ofReal
                  (s.norm_colTruncOf_sq_le N α l (hCC l))
            _ = ENNReal.ofReal (CC l) := by simp
        calc (∫⁻ α, ENNReal.ofReal (‖q α‖ ^ 2) ∂(gaussianMatrix r (d N)))
              * ENNReal.ofReal (1 / (d N : ℝ))
            ≤ ENNReal.ofReal (CC l) * ENNReal.ofReal (1 / (d N : ℝ)) :=
              mul_le_mul' hqb (le_refl _)
          _ = ENNReal.ofReal (CC l / (d N : ℝ)) := by
              rw [← ENNReal.ofReal_mul (hCC l), mul_one_div]
      calc ∑ l, ∫⁻ ξ, ENNReal.ofReal (overlapIdx ξ.2 a (s.colTruncOf N ξ.1 l (CC l)))
              ∂((gaussianMatrix r (d N)).prod (gaussianMatrix pN (d N)))
          ≤ ∑ l, ENNReal.ofReal (CC l / (d N : ℝ)) := Finset.sum_le_sum fun l _ => hone l
        _ = ENNReal.ofReal ((∑ l, CC l) / (d N : ℝ)) := by
            rw [Finset.sum_div, ENNReal.ofReal_sum_of_nonneg
              (fun l _ => div_nonneg (hCC l) hdpos.le)]
    calc ∑ a ∈ Finset.range r, ∫⁻ ξ,
            ∑ l, ENNReal.ofReal (overlapIdx ξ.2 a (s.colTruncOf N ξ.1 l (CC l)))
            ∂((gaussianMatrix r (d N)).prod (gaussianMatrix pN (d N)))
        ≤ ∑ _a ∈ Finset.range r, ENNReal.ofReal ((∑ l, CC l) / (d N : ℝ)) :=
          Finset.sum_le_sum hin
      _ = ENNReal.ofReal ((r : ℝ) * (∑ l, CC l) / (d N : ℝ)) := by
          rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul,
            ← ENNReal.ofReal_natCast r, ← ENNReal.ofReal_mul (Nat.cast_nonneg r),
            mul_div_assoc]
  rw [hint]
  refine (ENNReal.div_le_div_right hbound _).trans (le_of_eq ?_)
  rw [← ENNReal.ofReal_div_of_pos hκ, div_div]


end Deloc

/-! ### 6. The symmetry event (G7) -/

section Symmetry

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-- The spike directions are unit vectors. -/
theorem norm_spikeVec (s : RankRStack μ ns d r) (N : ℕ) (j : Fin r) :
    ‖s.spikeVec j N‖ = 1 := by
  have h := s.inner_spikeVec N j j
  rw [if_pos rfl] at h
  have h2 : ‖s.spikeVec j N‖ ^ 2 = 1 := by
    rw [← real_inner_self_eq_norm_sq]
    exact h
  calc ‖s.spikeVec j N‖ = Real.sqrt (‖s.spikeVec j N‖ ^ 2) :=
        (Real.sqrt_sq (norm_nonneg _)).symm
    _ = 1 := by rw [h2, Real.sqrt_one]

theorem norm_rvecOf_le (s : RankRStack μ ns d r) (N : ℕ) (j : Fin r)
    (α : Matrix (Fin r) (Fin (d N)) ℝ) : ‖s.rvecOf N j α‖ ≤ 1 := by
  by_cases h : IsUnit ((s.qmatOf N α)ᵀ * s.qmatOf N α).det
  · rw [s.rvecOf_of_isUnit N j h]
    have hle := DecompR.norm_toLp_rvecOfR_le (s.qmatOf N α)
      (v := WithLp.ofLp (s.spikeVec j N)) h
    rwa [WithLp.toLp_ofLp, s.norm_spikeVec N j] at hle
  · rw [s.rvecOf_of_not_isUnit N j h, norm_zero]
    norm_num

theorem dotProduct_qmatOf_rvecOf (s : RankRStack μ ns d r) (N : ℕ) (j : Fin r)
    (α : Matrix (Fin r) (Fin (d N)) ℝ) (k : Fin r) :
    ((s.qmatOf N α)ᵀ) k ⬝ᵥ WithLp.ofLp (s.rvecOf N j α) = 0 := by
  by_cases h : IsUnit ((s.qmatOf N α)ᵀ * s.qmatOf N α).det
  · rw [s.rvecOf_of_isUnit N j h, WithLp.ofLp_toLp]
    exact congrFun (DecompR.transpose_mulVec_rvecOfR (s.qmatOf N α)
      (v := WithLp.ofLp (s.spikeVec j N)) h) k
  · rw [s.rvecOf_of_not_isUnit N j h]
    simp

theorem measurable_transpose_qmatOf (s : RankRStack μ ns d r) (N : ℕ) :
    Measurable fun α : Matrix (Fin r) (Fin (d N)) ℝ => (s.qmatOf N α)ᵀ :=
  measurable_pi_lambda _ fun l => measurable_pi_lambda _ fun i =>
    s.measurable_qmatOf_apply N i l

/-- **(G7), Markov form.** The residual direction of the decomposition has small overlap with
the top-`r` eigenvectors of the Gram matrix. The bound is item Sym's `r / (d - r)` divided by
the level. -/
theorem measure_residual_ge_le (s : RankRStack μ ns d r) (hG : s.GaussianNoise) (N : ℕ)
    {U : Matrix (Fin (ns N)) (Fin r) ℝ} {pN : ℕ}
    {B : Ω N → Matrix (Fin pN) (Fin (d N)) ℝ}
    (hgram : ∀ ω, s.gram N ω = s.rankRW0 N ω U + s.qmatR N ω U * (s.qmatR N ω U)ᵀ)
    (hW0 : ∀ ω, s.rankRW0 N ω U = ((d N : ℝ))⁻¹ • ((B ω)ᵀ * B ω))
    (hlaw : HasLaw (fun ω => (Uᵀ * s.Zu N ω, B ω))
      ((gaussianMatrix r (d N)).prod (gaussianMatrix pN (d N))) (μ N))
    (hns0 : 0 < ns N) (hrk : r ≤ min (ns N) (d N)) (hrd : r < d N) (j : Fin r)
    {ξ0 : ℝ} (hξ : 0 < ξ0) :
    μ N {ω | ξ0 ≤ ∑ k ∈ Finset.range r,
        overlapIdx (s.X N ω) k (s.rvecOf N j (Uᵀ * s.Zu N ω))}
      ≤ ENNReal.ofReal ((r : ℝ) / (((d N : ℝ) - (r : ℝ)) * ξ0)) := by
  classical
  set G : Matrix (Fin r) (Fin (d N)) ℝ × Matrix (Fin pN) (Fin (d N)) ℝ → ℝ≥0∞ :=
    fun ζ => ∑ k ∈ Finset.range r, ENNReal.ofReal
      (overlapIdx (EdgeGlueR.ymatFlat (s.qmatOf N ζ.1)ᵀ (Real.sqrt (d N))⁻¹ ζ.2) k
        (s.rvecOf N j ζ.1)) with hGdef
  have hGm : Measurable G := by
    rw [hGdef]
    refine Finset.measurable_sum _ fun k _ => ?_
    exact EdgeGlueR.measurable_ofReal_overlapIdx_ymatFlat k (Real.sqrt (d N))⁻¹
      (Y := fun α => (s.qmatOf N α)ᵀ) (s.measurable_transpose_qmatOf N)
      (g := fun α => s.rvecOf N j α) (s.measurable_rvecOf N j)
  have hpt : ∀ ω, ENNReal.ofReal (∑ k ∈ Finset.range r,
      overlapIdx (s.X N ω) k (s.rvecOf N j (Uᵀ * s.Zu N ω)))
      = G (Uᵀ * s.Zu N ω, B ω) := by
    intro ω
    have hEq : (s.X N ω)ᵀ * s.X N ω
        = (EdgeGlueR.ymatFlat (s.qmatOf N (Uᵀ * s.Zu N ω))ᵀ (Real.sqrt (d N))⁻¹ (B ω))ᵀ
          * EdgeGlueR.ymatFlat (s.qmatOf N (Uᵀ * s.Zu N ω))ᵀ (Real.sqrt (d N))⁻¹ (B ω) := by
      rw [EdgeGlueR.gram_ymatFlat]
      exact s.gram_eq_YmatR N ω (hgram ω) (hW0 ω)
    rw [hGdef, ENNReal.ofReal_sum_of_nonneg
      (fun k _ => overlapIdx_nonneg _ k _)]
    refine Finset.sum_congr rfl fun k _ => ?_
    congr 1
    rw [overlapIdx, overlapIdx, specProjIdx_congr hEq (isHermitian_transpose_mul_self _)
      (isHermitian_transpose_mul_self _) k]
  have hmeas : AEMeasurable (fun ω => ENNReal.ofReal (∑ k ∈ Finset.range r,
      overlapIdx (s.X N ω) k (s.rvecOf N j (Uᵀ * s.Zu N ω)))) (μ N) := by
    refine AEMeasurable.congr ?_ (Filter.Eventually.of_forall fun ω => (hpt ω).symm)
    exact hGm.aemeasurable.comp_aemeasurable hlaw.aemeasurable
  have hset : {ω | ξ0 ≤ ∑ k ∈ Finset.range r,
        overlapIdx (s.X N ω) k (s.rvecOf N j (Uᵀ * s.Zu N ω))}
      = {ω | ENNReal.ofReal ξ0 ≤ ENNReal.ofReal (∑ k ∈ Finset.range r,
          overlapIdx (s.X N ω) k (s.rvecOf N j (Uᵀ * s.Zu N ω)))} := by
    ext ω
    simp only [Set.mem_ofPred_eq]
    rw [ENNReal.ofReal_le_ofReal_iff
      (Finset.sum_nonneg fun k _ => overlapIdx_nonneg _ k _)]
  have hξne : ENNReal.ofReal ξ0 ≠ 0 := by simp [ENNReal.ofReal_eq_zero, not_le, hξ]
  rw [hset]
  refine (meas_ge_le_lintegral_div hmeas hξne ENNReal.ofReal_ne_top).trans ?_
  have hint : ∫⁻ ω, ENNReal.ofReal (∑ k ∈ Finset.range r,
      overlapIdx (s.X N ω) k (s.rvecOf N j (Uᵀ * s.Zu N ω))) ∂(μ N)
      = ∫⁻ ζ, G ζ ∂((gaussianMatrix r (d N)).prod (gaussianMatrix pN (d N))) := by
    rw [lintegral_congr hpt]
    exact hlaw.lintegral_comp hGm.aemeasurable
  have hbound : ∫⁻ ζ, G ζ ∂((gaussianMatrix r (d N)).prod (gaussianMatrix pN (d N)))
      ≤ ENNReal.ofReal ((r : ℝ) / ((d N : ℝ) - (r : ℝ))) := by
    rw [lintegral_prod G hGm.aemeasurable]
    have hinner : ∀ᵐ α ∂(gaussianMatrix r (d N)),
        ∫⁻ Bm, G (α, Bm) ∂(gaussianMatrix pN (d N))
          ≤ ENNReal.ofReal ((r : ℝ) / ((d N : ℝ) - (r : ℝ))) := by
      filter_upwards [s.ae_ae_simpleSpec_YmatR hG N hgram hW0 hlaw hns0 hrk] with α hα
      have hcong : ∀ᵐ Bm ∂(gaussianMatrix pN (d N)), G (α, Bm)
          = ENNReal.ofReal (‖specProjTop
              ((YmatR (s.qmatOf N α)ᵀ (Real.sqrt (d N))⁻¹ Bm)ᵀ
                * YmatR (s.qmatOf N α)ᵀ (Real.sqrt (d N))⁻¹ Bm)
              (isHermitian_transpose_mul_self' _) r (s.rvecOf N j α)‖ ^ 2) := by
        filter_upwards [hα] with Bm hBm
        rw [hGdef, normSq_specProjTop_eq_sum_of_simpleSpec hBm (le_of_lt hrd),
          ENNReal.ofReal_sum_of_nonneg (fun k _ => sq_nonneg _)]
        refine Finset.sum_congr rfl fun k _ => ?_
        congr 1
        rw [overlapIdx, specProjIdx_congr (EdgeGlueR.gram_ymatFlat (s.qmatOf N α)ᵀ
          (Real.sqrt (d N))⁻¹ Bm) (isHermitian_transpose_mul_self _)
          (isHermitian_transpose_mul_self' _) k]
      rw [lintegral_congr_ae hcong]
      exact lintegral_normSq_specProjTop_YmatR_le hrd (s.qmatOf N α)ᵀ (Real.sqrt (d N))⁻¹ hα
        (s.norm_rvecOf_le N j α) (fun k => s.dotProduct_qmatOf_rvecOf N j α k)
    calc ∫⁻ α, ∫⁻ Bm, G (α, Bm) ∂(gaussianMatrix pN (d N)) ∂(gaussianMatrix r (d N))
        ≤ ∫⁻ _, ENNReal.ofReal ((r : ℝ) / ((d N : ℝ) - (r : ℝ)))
            ∂(gaussianMatrix r (d N)) := lintegral_mono_ae hinner
      _ = ENNReal.ofReal ((r : ℝ) / ((d N : ℝ) - (r : ℝ))) := by simp
  rw [hint]
  refine (ENNReal.div_le_div_right hbound _).trans (le_of_eq ?_)
  rw [← ENNReal.ofReal_div_of_pos hξ, div_div]


end Symmetry

/-! ### 7. The main theorem -/

section Main

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-- The scalar arithmetic of the assembly: with `m' > 64 r / (δ (-m₁))` the frame half of the
bound is at most `δ / 2`, and the residual half is `2 (δ / 4)`. -/
theorem edge_arith {rr : ℕ} (hrr : 0 < rr) {δ mder mm : ℝ} (hδ : 0 < δ) (hmm : 0 < mm)
    (hK : 64 * (rr : ℝ) / (δ * mm) < mder) :
    2 * ((rr : ℝ) / ((mder / 4) * (mm / 4))) + 2 * (δ / 4) ≤ δ := by
  have hrrR : (0 : ℝ) < rr := by exact_mod_cast hrr
  have hdm : 0 < δ * mm := mul_pos hδ hmm
  have h1 : 64 * (rr : ℝ) < mder * (δ * mm) := by
    rw [div_lt_iff₀ hdm] at hK
    linarith
  have hmder : 0 < mder := by nlinarith
  have hpos : (0 : ℝ) < (mder / 4) * (mm / 4) := by positivity
  have hkey : (rr : ℝ) / ((mder / 4) * (mm / 4)) ≤ δ / 4 := by
    rw [div_le_iff₀ hpos]
    nlinarith [h1]
  linarith

/-- **Task U7c, gap G1.** For every spike `j` and every `δ > 0` there is an `ε₀ > 0` such that
every `ε₁ ∈ (0, ε₀]` makes the projection of the spike direction on the top-`r` eigenvalues at
or below `bulkEdge c + ε₁` of squared norm more than `δ` with probability tending to `0`.

No supercriticality is assumed, so the statement holds for every spike of the mixed regime.
The consumer picks `ε₁ = min ε₀ ((min_k ρ_k - bulkEdge c) / 2)` and adds the count event
`tendsto_measure_eigenvalues₀_le` (`RankR/RMT/EdgeR.lean`) and the split
`Frame.norm_sq_specProjTop_split`. -/
theorem tendsto_measure_normSq_specProj_edge_gt [∀ N, IsProbabilityMeasure (μ N)]
    (s : RankRStack μ ns d r) (hG : s.GaussianNoise) {c : ℝ} (hc : 0 < c)
    (hdtop : Tendsto d atTop atTop)
    (hns : Tendsto (fun N => (ns N : ℝ) / d N) atTop (𝓝 c))
    {pp : ℕ → ℕ} (hn : ∀ N, ns N = r + pp N) (hp : ∀ N, 0 < pp N) (j : Fin r) :
    ∀ δ > 0, ∃ ε₀ > 0, ∀ ε₁, 0 < ε₁ → ε₁ ≤ ε₀ →
      Tendsto (fun N => μ N {ω |
        δ < ‖specProj (s.gram N ω)
          (topEigSet (s.gram N ω) (s.isHermitian_gram N ω) r ∩ Set.Iic (bulkEdge c + ε₁))
          (s.spikeVec j N)‖ ^ 2}) atTop (𝓝 0) := by
  classical
  intro δ hδ
  obtain ⟨U, hU, hgram, h, hBex⟩ := s.resolventLimitsR_of_gaussian_pair hG hc hdtop hns hn hp
  choose B hW0 hlaw using hBex
  -- 1. the scalars
  have hrpos : 0 < r := j.pos
  have hrR : (0 : ℝ) < r := by exact_mod_cast hrpos
  have hrne : (r : ℝ) ≠ 0 := ne_of_gt hrR
  have hb1 : (1 : ℝ) ≤ bulkEdge c := by
    have h0 : 0 ≤ Real.sqrt c := Real.sqrt_nonneg c
    rw [bulkEdge]
    nlinarith
  have hz₁ : bulkEdge c < bulkEdge c + 1 := by linarith
  have hm₁ : MP.m c (bulkEdge c + 1) < 0 := (MP.m_mem_Ioo hc hz₁).2
  have hmm : 0 < -MP.m c (bulkEdge c + 1) := by linarith
  obtain ⟨z₀, hz₀lo, hz₀hi, hz₀K⟩ :=
    R3minus.exists_mDeriv_gt hc (64 * (r : ℝ) / (δ * -MP.m c (bulkEdge c + 1)))
  have hmder : 0 < MP.mDeriv c z₀ := MP.mDeriv_pos hc hz₀lo
  set ε₀ : ℝ := (z₀ - bulkEdge c) / 2 with hε₀def
  have hε₀ : 0 < ε₀ := by rw [hε₀def]; linarith
  have hε₀ne : ε₀ ≠ 0 := ne_of_gt hε₀
  have hz₀eq : bulkEdge c + 2 * ε₀ = z₀ := by rw [hε₀def]; ring
  have hz₀pos : 0 < z₀ := by linarith
  set κ₀ : ℝ := ε₀ ^ 2 * MP.mDeriv c z₀ / (4 * r) with hκ₀def
  have hκ₀ : 0 < κ₀ := by rw [hκ₀def]; positivity
  set CC : Fin r → ℝ :=
    fun l => z₀ ^ 2 * ((s.coreEig l + 1) * MP.mDeriv c z₀ + MP.mDeriv c z₀ / (2 * r))
    with hCCdef
  have hCC0 : ∀ l, 0 ≤ CC l := by
    intro l
    have hl := s.coreEig_nonneg l
    simp only [hCCdef]
    positivity
  have hbreq : MP.mDeriv c z₀ / 2 - (r : ℝ) * κ₀ / ε₀ ^ 2 = MP.mDeriv c z₀ / 4 := by
    simp only [hκ₀def]
    field_simp
    ring
  refine ⟨ε₀, hε₀, ?_⟩
  intro ε₁ hε₁ hε₁le
  -- 2. the eventual size conditions
  have hr' : ∀ N, r ≤ ns N := fun N => by rw [hn N]; exact Nat.le_add_right r (pp N)
  have hppN : ∀ N, ns N - r = pp N := fun N => by rw [hn N]; omega
  have hcN : Tendsto (fun N => (pp N : ℝ) / d N) atTop (𝓝 c) := by
    have hh := tendsto_sub_ratio (nf := ns) r hr' hdtop hns
    simpa only [hppN] using hh
  have hdR : Tendsto (fun N => (d N : ℝ)) atTop atTop := tendsto_natCast_atTop_atTop.comp hdtop
  have hppinf : ∀ᶠ N in atTop, r ≤ pp N := by
    have h1 : ∀ᶠ N in atTop, c / 2 < (pp N : ℝ) / (d N : ℝ) :=
      hcN.eventually (eventually_gt_nhds (by linarith : c / 2 < c))
    have h2 : ∀ᶠ N in atTop, 2 * (r : ℝ) < c * (d N : ℝ) :=
      (Filter.Tendsto.const_mul_atTop hc hdR).eventually_gt_atTop (2 * (r : ℝ))
    filter_upwards [h1, h2] with N hn1 hn2
    have hdN : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast s.hd N
    rw [lt_div_iff₀ hdN] at hn1
    have hlt : (r : ℝ) < (pp N : ℝ) := by linarith
    exact_mod_cast hlt.le
  have hrdev : ∀ᶠ N in atTop, r < d N := hdtop.eventually_gt_atTop r
  -- 3. the bad families
  set bad1 : ∀ N, Set (Ω N) := fun N =>
    {ω | lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
      ≤ bulkEdge c + ε₀}ᶜ with hbad1def
  set bad1' : ∀ N, Set (Ω N) := fun N =>
    {ω | lamMax (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N))
      ≤ bulkEdge c + 1 / 2}ᶜ with hbad1'def
  set bad2 : ∀ N, Set (Ω N) := fun N =>
    {ω | SimpleSpec (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N)) r}ᶜ with hbad2def
  set bad3 : ∀ N, Set (Ω N) := fun N =>
    {ω | SimpleSpec (s.gram N ω) (s.isHermitian_gram N ω) r}ᶜ with hbad3def
  set bad4 : ∀ N, Set (Ω N) := fun N =>
    {ω | κ₀ ≤ ∑ a ∈ Finset.range r, ∑ l,
      ‖specProjIdx (s.rankRW0 N ω (U N)) (s.isHermitian_rankRW0 N ω (U N)) a
        (s.colTruncOf N ((U N)ᵀ * s.Zu N ω) l (CC l))‖ ^ 2} with hbad4def
  set bad5 : ∀ N, Set (Ω N) := fun N => ⋃ q : Fin r × Fin r,
    {ω | MP.mDeriv c z₀ / (2 * r) ≤ |R4.cform2 (s.rankRW0 N ω (U N)) z₀
      (fun i => s.qmatR N ω (U N) i q.1) (fun i => s.qmatR N ω (U N) i q.2)
      - (if q.1 = q.2 then (s.coreEig q.2 + 1) * MP.mDeriv c z₀ else 0)|} with hbad5def
  set bad6 : ∀ N, Set (Ω N) := fun N => ⋃ q : Fin r × Fin r,
    {ω | -MP.m c (bulkEdge c + 1) / (2 * r) ≤ |R4.cform (s.rankRW0 N ω (U N)) (bulkEdge c + 1)
      (fun i => s.qmatR N ω (U N) i q.1) (fun i => s.qmatR N ω (U N) i q.2)
      - (if q.1 = q.2 then (s.coreEig q.2 + 1) * MP.m c (bulkEdge c + 1) else 0)|}
    with hbad6def
  set bad7 : ∀ N, Set (Ω N) := fun N =>
    {ω | δ / 4 ≤ ∑ k ∈ Finset.range r,
      overlapIdx (s.X N ω) k (s.rvecOf N j ((U N)ᵀ * s.Zu N ω))} with hbad7def
  set bad : ∀ N, Set (Ω N) := fun N =>
    bad1 N ∪ (bad1' N ∪ (bad2 N ∪ (bad3 N ∪ (bad4 N ∪ (bad5 N ∪ (bad6 N ∪ bad7 N))))))
    with hbaddef
  -- 4. each bad family vanishes
  have hv1 : Tendsto (fun N => μ N (bad1 N)) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (s.measurableSet_lamMax_rankRW0_le N (U N) _).nullMeasurableSet)
      (h.edge ε₀ hε₀)
  have hv1' : Tendsto (fun N => μ N (bad1' N)) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (s.measurableSet_lamMax_rankRW0_le N (U N) _).nullMeasurableSet)
      (h.edge (1 / 2) (by norm_num))
  have hv2 : Tendsto (fun N => μ N (bad2 N)) atTop (𝓝 0) := by
    refine EdgeGlueR.tendsto_measure_zero_of_eventually_zero ?_
    filter_upwards [hppinf, hrdev] with N hpN hrN
    exact ae_iff.mp (s.ae_simpleSpec_rankRW0 N (hW0 N)
      (measurePreserving_snd.fun_comp_hasLaw (hlaw N)) (hp N) (le_min hpN hrN.le))
  have hv3 : Tendsto (fun N => μ N (bad3 N)) atTop (𝓝 0) := by
    refine EdgeGlueR.tendsto_measure_zero_of_eventually_zero ?_
    filter_upwards [hrdev] with N hrN
    have hns0 : 0 < ns N := by rw [hn N]; omega
    exact ae_iff.mp (s.ae_simpleSpec_gram hG N hns0 (le_min (hr' N) hrN.le))
  have hv4 : Tendsto (fun N => μ N (bad4 N)) atTop (𝓝 0) := by
    have hreal : Tendsto (fun N => (r : ℝ) * (∑ l, CC l) / ((d N : ℝ) * κ₀)) atTop (𝓝 0) := by
      have h1 : Tendsto (fun N => (d N : ℝ) * κ₀) atTop atTop :=
        Filter.Tendsto.atTop_mul_const hκ₀ hdR
      have h3 := h1.inv_tendsto_atTop.const_mul ((r : ℝ) * (∑ l, CC l))
      simpa [div_eq_mul_inv] using h3
    have htend : Tendsto (fun N => ENNReal.ofReal
        ((r : ℝ) * (∑ l, CC l) / ((d N : ℝ) * κ₀))) atTop (𝓝 0) := by
      have h3 := (ENNReal.continuous_ofReal.tendsto 0).comp hreal
      rw [ENNReal.ofReal_zero] at h3
      exact h3
    refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds htend
      (Filter.Eventually.of_forall fun _ => zero_le) ?_
    filter_upwards [hppinf, hrdev] with N hpN hrN
    exact s.measure_kappa_ge_le N (hW0 N) (hlaw N) (le_min hpN hrN.le) hCC0 hκ₀
  have hv5 : Tendsto (fun N => μ N (bad5 N)) atTop (𝓝 0) :=
    tendsto_measure_zero_iUnion fun q =>
      s.tendstoInProb_cform2_qmatR h q.1 q.2 hz₀lo (MP.mDeriv c z₀ / (2 * r)) (by positivity)
  have hv6 : Tendsto (fun N => μ N (bad6 N)) atTop (𝓝 0) :=
    tendsto_measure_zero_iUnion fun q =>
      s.tendstoInProb_cform_qmatR h q.1 q.2 hz₁
        (-MP.m c (bulkEdge c + 1) / (2 * r)) (by positivity)
  have hv7 : Tendsto (fun N => μ N (bad7 N)) atTop (𝓝 0) := by
    have hreal : Tendsto (fun N => (r : ℝ) / (((d N : ℝ) - (r : ℝ)) * (δ / 4)))
        atTop (𝓝 0) := by
      have h1 : Tendsto (fun N => ((d N : ℝ) - (r : ℝ)) * (δ / 4)) atTop atTop :=
        Filter.Tendsto.atTop_mul_const (by linarith : (0 : ℝ) < δ / 4)
          ((tendsto_natCast_atTop_atTop.comp hdtop).atTop_add tendsto_const_nhds)
      have h3 := h1.inv_tendsto_atTop.const_mul (r : ℝ)
      simpa [div_eq_mul_inv] using h3
    have htend : Tendsto (fun N => ENNReal.ofReal
        ((r : ℝ) / (((d N : ℝ) - (r : ℝ)) * (δ / 4)))) atTop (𝓝 0) := by
      have h3 := (ENNReal.continuous_ofReal.tendsto 0).comp hreal
      rw [ENNReal.ofReal_zero] at h3
      exact h3
    refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds htend
      (Filter.Eventually.of_forall fun _ => zero_le) ?_
    filter_upwards [hrdev] with N hrN
    have hns0 : 0 < ns N := by rw [hn N]; omega
    exact s.measure_residual_ge_le hG N (hgram N) (hW0 N) (hlaw N) hns0
      (le_min (hr' N) hrN.le) hrN j (by linarith)
  have hvbad : Tendsto (fun N => μ N (bad N)) atTop (𝓝 0) :=
    tendsto_measure_zero_union hv1 (tendsto_measure_zero_union hv1'
      (tendsto_measure_zero_union hv2 (tendsto_measure_zero_union hv3
        (tendsto_measure_zero_union hv4 (tendsto_measure_zero_union hv5
          (tendsto_measure_zero_union hv6 hv7))))))
  -- 5. the deterministic step on the intersection of the good events
  refine OutliersR.tendsto_measure_zero_of_eventually_subset (t := bad) ?_ hvbad
  filter_upwards [hppinf, hrdev] with N hpN hrN
  intro ω hω
  by_contra hbad
  simp only [hbaddef, hbad1def, hbad1'def, hbad2def, hbad3def, hbad4def, hbad5def, hbad6def,
    hbad7def, Set.mem_union, Set.mem_compl_iff, Set.mem_iUnion, Set.mem_ofPred_eq,
    not_or, not_exists, not_le, not_lt, not_not] at hbad
  obtain ⟨he1, he1', he2, he3, he4, he5, he6, he7⟩ := hbad
  have hWh := s.isHermitian_rankRW0 N ω (U N)
  have hSh := s.isHermitian_gram N ω
  have hpsd := s.eigenvalues_rankRW0_nonneg N ω (U N)
  have hqeq : s.qmatOf N ((U N)ᵀ * s.Zu N ω) = s.qmatR N ω (U N) := s.qmatOf_eq N ω (U N)
  have hlamlt0 : lamMax (s.rankRW0 N ω (U N)) hWh < z₀ := by
    have hh : bulkEdge c + ε₀ < z₀ := by linarith
    linarith
  have hlamlt1 : lamMax (s.rankRW0 N ω (U N)) hWh < bulkEdge c + 1 := by linarith
  -- (i) the `μmin` bound from G5
  have hclose5 : ∀ k l : Fin r,
      |R4.cform2 (s.rankRW0 N ω (U N)) z₀ (fun i => s.qmatR N ω (U N) i k)
          (fun i => s.qmatR N ω (U N) i l)
        - (if k = l then (s.coreEig k + 1) * MP.mDeriv c z₀ else 0)|
        ≤ MP.mDeriv c z₀ / (2 * r) := by
    intro k l
    rcases eq_or_ne k l with rfl | hne
    · exact (he5 (k, k)).le
    · rw [if_neg hne]
      have hq := he5 (k, l)
      rw [if_neg hne] at hq
      exact hq.le
  have hmin : ∀ y : Fin r → ℝ, (MP.mDeriv c z₀ / 2) * (y ⬝ᵥ y)
      ≤ R4.qform2 (s.rankRW0 N ω (U N)) z₀ (s.qmatR N ω (U N) *ᵥ y) := by
    intro y
    rw [EdgeGlueDetR.qform2_mulVec_expand]
    exact EdgeGlueDetR.le_sum_of_close_to_diag hrpos s.coreEig_nonneg hmder.le hclose5 y
  -- (ii) the column norms
  have hcolnorm : ∀ l : Fin r,
      (fun i => s.qmatR N ω (U N) i l) ⬝ᵥ (fun i => s.qmatR N ω (U N) i l) ≤ CC l := by
    intro l
    have h1 := EdgeGlueDetR.le_qform2_of_psd hWh hpsd hz₀pos hlamlt0
      (fun i => s.qmatR N ω (U N) i l)
    have hq := he5 (l, l)
    rw [if_pos rfl] at hq
    have hq2 := (abs_lt.mp hq).2
    rw [R4.cform2_self] at hq2
    have hz2 : (0 : ℝ) < z₀ ^ 2 := pow_pos hz₀pos 2
    rw [div_le_iff₀ hz2] at h1
    simp only [hCCdef]
    nlinarith [h1, hq2, hz2]
  have hcolT : ∀ l : Fin r, s.colTruncOf N ((U N)ᵀ * s.Zu N ω) l (CC l)
      = WithLp.toLp 2 (fun i => s.qmatR N ω (U N) i l) := by
    intro l
    have hcol : s.colOf N ((U N)ᵀ * s.Zu N ω) l
        = WithLp.toLp 2 (fun i => s.qmatR N ω (U N) i l) := by
      rw [RankRStack.colOf, hqeq]
    have hnorm : ‖s.colOf N ((U N)ᵀ * s.Zu N ω) l‖ ^ 2 ≤ CC l := by
      rw [hcol, EdgeDetR.norm_toLp_sq]
      exact hcolnorm l
    rw [s.colTruncOf_eq_colOf N _ l hnorm, hcol]
  have hκ : ∑ a ∈ Finset.range r, ∑ l,
      ‖specProjIdx (s.rankRW0 N ω (U N)) hWh a
        (WithLp.toLp 2 fun i => s.qmatR N ω (U N) i l)‖ ^ 2 ≤ κ₀ := by
    have hh := he4
    simp only [hcolT] at hh
    exact hh.le
  -- (iii) the `μQ` bound from G6
  have hclose6 : ∀ k l : Fin r,
      |(-R4.cform (s.rankRW0 N ω (U N)) (bulkEdge c + 1) (fun i => s.qmatR N ω (U N) i k)
          (fun i => s.qmatR N ω (U N) i l))
        - (if k = l then (s.coreEig k + 1) * -MP.m c (bulkEdge c + 1) else 0)|
        ≤ -MP.m c (bulkEdge c + 1) / (2 * r) := by
    intro k l
    rcases eq_or_ne k l with rfl | hne
    · have hq := he6 (k, k)
      rw [if_pos rfl] at hq
      rw [if_pos rfl, show (-R4.cform (s.rankRW0 N ω (U N)) (bulkEdge c + 1)
            (fun i => s.qmatR N ω (U N) i k) (fun i => s.qmatR N ω (U N) i k))
          - (s.coreEig k + 1) * -MP.m c (bulkEdge c + 1)
          = -(R4.cform (s.rankRW0 N ω (U N)) (bulkEdge c + 1)
            (fun i => s.qmatR N ω (U N) i k) (fun i => s.qmatR N ω (U N) i k)
            - (s.coreEig k + 1) * MP.m c (bulkEdge c + 1)) from by ring, abs_neg]
      exact hq.le
    · have hq := he6 (k, l)
      rw [if_neg hne] at hq
      rw [if_neg hne, show (-R4.cform (s.rankRW0 N ω (U N)) (bulkEdge c + 1)
            (fun i => s.qmatR N ω (U N) i k) (fun i => s.qmatR N ω (U N) i l)) - 0
          = -(R4.cform (s.rankRW0 N ω (U N)) (bulkEdge c + 1)
            (fun i => s.qmatR N ω (U N) i k) (fun i => s.qmatR N ω (U N) i l) - 0) from by
          ring, abs_neg]
      exact hq.le
  have hdom6 : ∀ y : Fin r → ℝ, (-MP.m c (bulkEdge c + 1) / 2) * (y ⬝ᵥ y)
      ≤ -R4.qform (s.rankRW0 N ω (U N)) (bulkEdge c + 1) (s.qmatR N ω (U N) *ᵥ y) := by
    intro y
    have hd := EdgeGlueDetR.le_sum_of_close_to_diag hrpos s.coreEig_nonneg hmm.le hclose6 y
    have hneg : ∑ k, ∑ l, y k * y l *
        (-R4.cform (s.rankRW0 N ω (U N)) (bulkEdge c + 1)
          (fun i => s.qmatR N ω (U N) i k) (fun i => s.qmatR N ω (U N) i l))
        = -R4.qform (s.rankRW0 N ω (U N)) (bulkEdge c + 1) (s.qmatR N ω (U N) *ᵥ y) := by
      rw [EdgeGlueDetR.qform_mulVec_expand, ← Finset.sum_neg_distrib]
      refine Finset.sum_congr rfl fun k _ => ?_
      rw [← Finset.sum_neg_distrib]
      exact Finset.sum_congr rfl fun l _ => by ring
    rw [← hneg]
    exact hd
  have hQQ : ∀ y : Fin r → ℝ, (-MP.m c (bulkEdge c + 1) / 4) * (y ⬝ᵥ y)
      ≤ y ⬝ᵥ (((s.qmatR N ω (U N))ᵀ * s.qmatR N ω (U N)) *ᵥ y) := by
    intro y
    rw [← EdgeGlueR.dotProduct_mulVec_self]
    have h1 := EdgeGlueDetR.mul_neg_qform_le hWh hlamlt1 (s.qmatR N ω (U N) *ᵥ y)
    have h2 := hdom6 y
    have hyy : 0 ≤ y ⬝ᵥ y := DecompR.dotProduct_self_nonneg' y
    have h3 : (1 : ℝ) / 2 ≤ bulkEdge c + 1 - lamMax (s.rankRW0 N ω (U N)) hWh := by linarith
    have h4 : 0 ≤ -R4.qform (s.rankRW0 N ω (U N)) (bulkEdge c + 1)
        (s.qmatR N ω (U N) *ᵥ y) := le_trans (by positivity) h2
    nlinarith [h1, h2, h3, h4, hyy]
  have hμQ : (0 : ℝ) < -MP.m c (bulkEdge c + 1) / 4 := by linarith
  have hunit : IsUnit ((s.qmatR N ω (U N))ᵀ * s.qmatR N ω (U N)).det :=
    EdgeGlueDetR.isUnit_det_of_lower_bound (s.qmatR N ω (U N)) hμQ hQQ
  -- (iv) the residual from G7
  have hrvec : s.rvecOf N j ((U N)ᵀ * s.Zu N ω)
      = WithLp.toLp 2 (DecompR.rvecOfR (s.qmatR N ω (U N))
          (WithLp.ofLp (s.spikeVec j N))) := by
    have hu : IsUnit ((s.qmatOf N ((U N)ᵀ * s.Zu N ω))ᵀ
        * s.qmatOf N ((U N)ᵀ * s.Zu N ω)).det := by rw [hqeq]; exact hunit
    rw [s.rvecOf_of_isUnit N j hu, hqeq]
  have hrest : ∑ k ∈ Finset.range r,
      ‖specProjIdx (s.gram N ω) hSh k (WithLp.toLp 2 (DecompR.rvecOfR (s.qmatR N ω (U N))
        (WithLp.ofLp (s.spikeVec j N))))‖ ^ 2 ≤ δ / 4 := by
    have hh := he7
    rw [hrvec] at hh
    exact hh.le
  -- (v) the assembly
  have hedge : lamMax (s.rankRW0 N ω (U N)) hWh + ε₀ ≤ z₀ := by linarith
  have hτ : bulkEdge c + ε₁ ≤ z₀ := by linarith
  have hdet := EdgeGlueDetR.normSq_specProj_edge_le hWh hSh (hgram N ω) he2 he3 hrN.le
    (z₀ := z₀) (ε₀ := ε₀) (τ := bulkEdge c + ε₁) (κ := κ₀)
    (μmin := MP.mDeriv c z₀ / 2) (μQ := -MP.m c (bulkEdge c + 1) / 4) (ξ := δ / 4)
    hε₀ hedge hτ hκ hmin (by rw [hbreq]; linarith) hμQ hQQ (s.norm_spikeVec N j) hrest
  rw [hbreq] at hdet
  have hfinal := hdet.trans (edge_arith hrpos hδ hmm hz₀K)
  exact absurd hω (not_lt.mpr hfinal)

/-- The good-set form of `tendsto_measure_normSq_specProj_edge_gt`, as the task states it. No
measurability of the good set is needed: it is the complement of the bad set. -/
theorem tendsto_measure_normSq_specProj_edge_le [∀ N, IsProbabilityMeasure (μ N)]
    (s : RankRStack μ ns d r) (hG : s.GaussianNoise) {c : ℝ} (hc : 0 < c)
    (hdtop : Tendsto d atTop atTop)
    (hns : Tendsto (fun N => (ns N : ℝ) / d N) atTop (𝓝 c))
    {pp : ℕ → ℕ} (hn : ∀ N, ns N = r + pp N) (hp : ∀ N, 0 < pp N) (j : Fin r) :
    ∀ δ > 0, ∃ ε₀ > 0, ∀ ε₁, 0 < ε₁ → ε₁ ≤ ε₀ →
      Tendsto (fun N => μ N {ω |
        ‖specProj (s.gram N ω)
          (topEigSet (s.gram N ω) (s.isHermitian_gram N ω) r ∩ Set.Iic (bulkEdge c + ε₁))
          (s.spikeVec j N)‖ ^ 2 ≤ δ}) atTop (𝓝 1) := by
  intro δ hδ
  obtain ⟨ε₀, hε₀, hmain⟩ :=
    s.tendsto_measure_normSq_specProj_edge_gt hG hc hdtop hns hn hp j δ hδ
  refine ⟨ε₀, hε₀, fun ε₁ hε₁ hε₁le => ?_⟩
  refine tendsto_measure_one_of_bad ?_ (hmain ε₁ hε₁ hε₁le)
  intro N ω hω
  simp only [Set.mem_compl_iff, Set.mem_ofPred_eq, not_le] at hω
  exact hω

end Main

end RankRStack

end StackedSVD
