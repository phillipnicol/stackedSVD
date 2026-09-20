/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.AlignedComponent
import StackedSVD.RankR.GeneralFrob
import StackedSVD.LinAlg.SpecIdxPerturb

/-!
# The component clause of `thm:rank_r_svdstack`

Track D item D4, steps 4 to 6 of `notes/archive/rankr_D4_plan.md`. Paper statement
`main_paper.tex:2405`: under the exactly aligned rank-`r` model,
`(v_jᵀ v̂_{j,svdstack})² → S_j / (S_j + 1)` for every component `j`, with
`S_j = ∑_i β_ij² / (1 - β_ij²)`. The companion aggregate clause is item D3
(`UnalignedModelR.thm_rank_r_svdstack_aggregate`, `RankR/WeightedUpperG.lean`).

Write `W⋆ = optWG β`, `A⋆ = W⋆ A_{β,R} W⋆ᵀ`, `C = W⋆ B_R`, `y_j` for column `j` of `Ṽ_{W⋆} V`
and `c_j` for column `j` of `C`.

## Content

1. **The frame-free quantity.** `compRGW W j = ‖P_j y_j‖² / λ_j(Ṽ_W Ṽ_Wᵀ)`, with `P_j` the
   projector at the sorted eigenvalue index `j`. It needs no eigenframe and no simplicity, so
   it is defined at every `ω`. Rank-one mirror: `overlap` (`Defs.lean`).
2. **The two branches.** At `0 < S_j` the limit matrix `A⋆` has gaps at `j` and at `j + 1`
   (`RankR/AlignedComponent.lean`), the quotient functional `compFun`
   (`LinAlg/SpecIdxPerturb.lean`) is continuous there, and the two entrywise limits `gramRWG`
   and `VtVRWG` (`RankR/GeneralFrob.lean`) transport it. At `S_j = 0` there is no gap and no
   simplicity: `c_j = 0`, so `‖y_j‖² → 0` while `λ_j → 1`, and a projector contracts, which
   squeezes the quotient to `0 = S_j / (S_j + 1)`.
3. **The paper-facing forms.** `thm_rank_r_svdstack_component_eig_of_sep` reads the estimator
   at the canonical frame `(topEigMat, topEigVal)` of `Ṽ_{W⋆} Ṽ_{W⋆}ᵀ`, which is sorted by
   construction, so it carries no frame hypothesis.
   `thm_rank_r_svdstack_component_frame_of_sep` takes any measurable frame that is sorted with
   probability tending to one. `IsTopEigFrame` alone does not order the columns, so the sorted
   conjunct is needed: a permuted frame permutes the conclusion (modeling choice M3).
4. **The hypothesis on `S`.** Every core takes `hsep : SaggSep β`, the separation the aligned
   spectrum reads. The paper's hypothesis (b) `S_j ≠ S_k` gives it
   (`Sagg_sep_of_strictAnti`), and the model gives it with no hypothesis at all
   (`saggSep_of_model`), which is `thm_rank_r_svdstack_component_of_model`. There is no
   hypothesis `0 < S_j` anywhere.

STATUS 2026-09-02: proved, 0 `sorry`.
-/

open MeasureTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

/-! ### 1. Deterministic bridges -/

section Bridges

variable {p r : ℕ} {A : Matrix (Fin p) (Fin p) ℝ}

/-- The squared Euclidean norm of a coordinate function. -/
theorem normSq_toLp {n : ℕ} (f : Fin n → ℝ) :
    ‖(WithLp.toLp 2 f : EuclideanSpace ℝ (Fin n))‖ ^ 2 = ∑ i, f i ^ 2 := by
  rw [← real_inner_self_eq_norm_sq, real_inner_eq_dotProduct, dotProduct]
  exact Finset.sum_congr rfl fun i _ => (sq (f i)).symm

/-- Column `j` of the canonical frame is the sorted eigenvector at index `j`. -/
theorem frameCol_topEigMat_eq_vEig (hA : A.IsHermitian) (hrp : r ≤ Fintype.card (Fin p))
    (j : Fin r) : frameCol (topEigMat hA hrp) j = vEig A hA (j : ℕ) := by
  rw [frameCol_topEigMat, vEig_eq A hA (lt_of_lt_of_le j.isLt hrp)]
  rfl

/-- The canonical frame eigenvalue at column `j` is `eigVal` at index `j`. -/
theorem topEigVal_eq_eigVal (hA : A.IsHermitian) (hrp : r ≤ Fintype.card (Fin p)) (j : Fin r) :
    topEigVal hA hrp j = eigVal A hA (j : ℕ) := by
  rw [topEigVal, eigVal_eq A hA (lt_of_lt_of_le j.isLt hrp)]
  rfl

/-- A unit eigenvector at the sorted index `j` computes the index projector, when that index
is simple. The sign is absorbed by the square. This is the general-frame companion of
`normSq_specProjIdx_eq_inner_sq` (`LinAlg/SpecIdxPerturb.lean`), which reads `vEig`. -/
theorem inner_sq_eq_normSq_specProjIdx {hA : A.IsHermitian} {j : ℕ}
    (hsimple : SimpleIdx A hA j) (hj : j < Fintype.card (Fin p))
    {q : EuclideanSpace ℝ (Fin p)} (hq : toOp A q = hA.eigenvalues₀ ⟨j, hj⟩ • q)
    (hqn : ‖q‖ = 1) (y : EuclideanSpace ℝ (Fin p)) :
    ⟪q, y⟫_ℝ ^ 2 = ‖specProjIdx A hA j y‖ ^ 2 := by
  have hmem : q ∈ specSpace A (eigSetIdx A hA j) := by
    rw [eigSetIdx_eq_singleton A hA hj, specSpace_singleton]
    exact Module.End.mem_eigenspace_iff.mpr hq
  have hproj : specProjIdx A hA j q = q :=
    Submodule.starProjection_eq_self_iff.mpr hmem
  have hvn : ‖vEig A hA j‖ = 1 := norm_vEig A hA hj
  set t := ⟪vEig A hA j, q⟫_ℝ with ht
  have ht2 : t ^ 2 = 1 := by
    have h := normSq_specProjIdx_eq_inner_sq hsimple q
    rw [hproj, hqn, one_pow] at h
    exact h.symm
  have hqv : ⟪q, vEig A hA j⟫_ℝ = t := by rw [ht, real_inner_comm]
  have hzero : ‖q - t • vEig A hA j‖ ^ 2 = 0 := by
    rw [norm_sub_sq_real, real_inner_smul_right, hqv, norm_smul, Real.norm_eq_abs, hvn,
      hqn, mul_one]
    nlinarith [sq_abs t]
  have hqe : q = t • vEig A hA j := by
    have := pow_eq_zero_iff (n := 2) (by norm_num) |>.mp hzero
    rw [norm_eq_zero, sub_eq_zero] at this
    exact this
  rw [hqe, real_inner_smul_left, mul_pow, ht2, one_mul,
    normSq_specProjIdx_eq_inner_sq hsimple y]

/-- A frame column of an orthonormal frame is a unit vector. -/
theorem norm_frameCol_of_ortho {Y : Matrix (Fin p) (Fin r) ℝ} (h : Yᵀ * Y = 1) (j : Fin r) :
    ‖frameCol Y j‖ = 1 := by
  have h2 : ‖frameCol Y j‖ ^ 2 = 1 := by
    rw [← real_inner_self_eq_norm_sq, inner_frameCol_of_ortho h, if_pos rfl]
  rw [← Real.sqrt_sq (norm_nonneg (frameCol Y j)), h2, Real.sqrt_one]

/-- Cauchy-Schwarz at a unit vector, divided by a nonnegative scalar. At `L = 0` both sides
are `0`, so no positivity of `L` is needed. -/
theorem inner_sq_div_le (v y : EuclideanSpace ℝ (Fin p)) (hvn : ‖v‖ = 1) {L : ℝ} (hL : 0 ≤ L) :
    |⟪v, y⟫_ℝ ^ 2 / L| ≤ ‖y‖ ^ 2 / L := by
  have habs : |⟪v, y⟫_ℝ| ≤ ‖y‖ := by
    have h := abs_real_inner_le_norm v y
    rwa [hvn, one_mul] at h
  have hnn : 0 ≤ ‖y‖ ^ 2 - ⟪v, y⟫_ℝ ^ 2 := by
    nlinarith [habs, sq_abs (⟪v, y⟫_ℝ), abs_nonneg (⟪v, y⟫_ℝ), norm_nonneg y]
  have hkey : 0 ≤ (‖y‖ ^ 2 - ⟪v, y⟫_ℝ ^ 2) / L := div_nonneg hnn hL
  rw [sub_div] at hkey
  rw [abs_of_nonneg (div_nonneg (sq_nonneg _) hL)]
  linarith

end Bridges

/-! ### 2. The limit eigenvalue at one component index -/

section AlignedEigVal

variable {M r : ℕ}

/-- `λ_j(A⋆) = 1 + S_j` at a component index, in the `eigVal` form. Item 2.4 of
`notes/archive/rankr_D4_plan.md` read at the free index `j`. -/
theorem eigVal_ABlockW_aligned (β : Fin M → Fin r → ℝ) (h0 : ∀ i j, 0 ≤ β i j)
    (h1 : ∀ i j, β i j < 1) (hanti : Antitone (Sagg β)) (hM : 0 < M) (j : Fin r) :
    eigVal (ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r))
        (isHermitian_ABlockW (optWG (rk := alignedRk M r) β) β (alignedR M r)) (j : ℕ)
      = 1 + Sagg β j := by
  have hjc : (j : ℕ) < Fintype.card (Fin (rtot (alignedRk M r))) :=
    lt_card_rtot_alignedRk hM j.isLt
  rw [eigVal_eq _ _ hjc, eigenvalues₀_ABlockW_optWG_aligned β h0 h1 hanti hM, dif_pos j.isLt,
    Fin.eta]

end AlignedEigVal

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-! ### 3. The two new objects -/

/-- Column `k` of `Ṽ_W V`, as a vector of `EuclideanSpace`. -/
noncomputable def VtVWGcol (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (N : ℕ) (ω : Ω N) (k : Fin r) :
    EuclideanSpace ℝ (Fin (rtot rk)) :=
  WithLp.toLp 2 fun p => m.VtVWG W N ω p k

/-- The paper's `(v_jᵀ v̂_{j,svdstack})²` at the weight matrix `W`, frame free. -/
noncomputable def compRGW (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (j : Fin r) (N : ℕ) (ω : Ω N) : ℝ :=
  ‖specProjIdx (m.gramWG W N ω) (m.isHermitian_gramWG W N ω) (j : ℕ)
      (m.VtVWGcol W N ω j)‖ ^ 2
    / eigVal (m.gramWG W N ω) (m.isHermitian_gramWG W N ω) (j : ℕ)

theorem frameCol_VtVWG (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (N : ℕ) (ω : Ω N) (k : Fin r) :
    frameCol (m.VtVWG W N ω) k = m.VtVWGcol W N ω k := rfl

/-- **The paper's frame quantity, squared** (`main_paper.tex:2405`). The diagonal entry `j` of
`V̂_svdstack(W)ᵀ V` is `(√λ_j)^{-1} ⟪q_j, y_j⟫` with `q_j` column `j` of the frame and `y_j`
column `j` of `Ṽ_W V`, so its square is `⟪q_j, y_j⟫² / λ_j`. The identity also holds at
`λ_j = 0`, where both sides are `0`. -/
theorem vhatSvdstackGW_entry_sq (m : UnalignedModelR μ M n d r rk)
    (W : Matrix (Fin (rtot rk)) (Fin (rtot rk)) ℝ) (N : ℕ) (ω : Ω N)
    (Q : Matrix (Fin (rtot rk)) (Fin r) ℝ) (lam : Fin r → ℝ) (j : Fin r)
    (hlam : 0 ≤ lam j) :
    (((m.vhatSvdstackGW W N ω Q lam)ᵀ * m.V N) j j) ^ 2
      = ⟪frameCol Q j, m.VtVWGcol W N ω j⟫_ℝ ^ 2 / lam j := by
  have hinv : (Real.sqrt (lam j))⁻¹ ^ 2 = (lam j)⁻¹ := by
    rw [inv_pow, Real.sq_sqrt hlam]
  rw [m.vhatSvdstackGW_transpose_mul_V W N ω Q lam, Matrix.diagonal_mul,
    ← m.frameCol_VtVWG W N ω j, inner_frameCol, mul_pow, hinv, inv_mul_eq_div]

end UnalignedModelR

/-! ### 4. The component clause -/

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-- `gramRWG` in the exactly aligned model: `Ṽ_{W⋆} Ṽ_{W⋆}ᵀ → A⋆` entrywise. -/
theorem gramWG_tendsto_aligned (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hR : ∀ i, m.R i = 1)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise)
    (p q : Fin (rtot (alignedRk M r))) :
    TendstoInProb μ (fun N ω => m.gramWG (optWG β) N ω p q)
      (ABlockW (optWG β) β (alignedR M r) p q) := by
  have hRe : m.R = alignedR M r := funext hR
  have h := m.gramRWG c β (optWG β) hβdef law hG p q
  rwa [hRe] at h

/-- `VtVRWG` in the exactly aligned model: `Ṽ_{W⋆} V → C` entrywise. -/
theorem VtVWG_tendsto_aligned (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hR : ∀ i, m.R i = 1)
    (law : ∀ i, (m.tbl i).TableLawR (c i))
    (p : Fin (rtot (alignedRk M r))) (k : Fin r) :
    TendstoInProb μ (fun N ω => m.VtVWG (optWG β) N ω p k)
      (BBlockW (optWG β) β (alignedR M r) p k) := by
  have hRe : m.R = alignedR M r := funext hR
  have h := m.VtVRWG c β (optWG β) hβdef law p k
  rwa [hRe] at h

/-- **`thm:rank_r_svdstack`, component clause, frame free** (`main_paper.tex:2405`, Track D
item D4). `compRGW` is `‖P_j y_j‖² / λ_j`, the paper's `(v_jᵀ v̂_j)²` written with no choice of
eigenframe; `thm_rank_r_svdstack_component_eig_of_sep` reads it at the canonical frame.

`hsep : SaggSep β` is the separation the aligned spectrum needs; the paper's hypothesis (b)
gives it through `Sagg_sep_of_strictAnti`, and the model gives it with no hypothesis at all
through `saggSep_of_model`. There is no hypothesis `0 < S_j`: at `S_j = 0` the numerator
vanishes and the limit is `0`, which is `S_j / (S_j + 1)`. -/
theorem thm_rank_r_svdstack_component_of_sep
    (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hsep : SaggSep β)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) (j : Fin r) :
    TendstoInProb μ (fun N ω => m.compRGW (optWG β) j N ω) (Sagg β j / (Sagg β j + 1)) := by
  have hβ01 : ∀ i k, 0 ≤ β i k ∧ β i k < 1 := fun i k => by
    rw [hβdef i k]
    exact beta_mem_Ico (hc i)
  have h0 : ∀ i k, 0 ≤ β i k := fun i k => (hβ01 i k).1
  have h1 : ∀ i k, β i k < 1 := fun i k => (hβ01 i k).2
  have hsq : ∀ i k, β i k ^ 2 < 1 := fun i k => by nlinarith [h0 i k, h1 i k]
  have hev := eigVal_ABlockW_aligned β h0 h1 hsep.1 hM j
  have hgram := m.gramWG_tendsto_aligned c β hβdef hR law hG
  have hVtV := m.VtVWG_tendsto_aligned c β hβdef hR law
  rcases eq_or_lt_of_le (Sagg_nonneg hsq j) with hzero | hpos
  · -- branch Z: `S_j = 0`, a squeeze; the limit matrix has no gap at `j`
    have hSz : Sagg β j = 0 := hzero.symm
    have hlimval : Sagg β j / (Sagg β j + 1) = 0 := by rw [hSz]; norm_num
    rw [hlimval]
    have hden : TendstoInProb μ
        (fun N ω => eigVal (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω)
          (j : ℕ)) 1 := by
      have hpi : TendstoInProbPi μ
          (fun N ω => fun t : Fin (rtot (alignedRk M r)) × Fin (rtot (alignedRk M r)) =>
            m.gramWG (optWG β) N ω t.1 t.2)
          (fun t : Fin (rtot (alignedRk M r)) × Fin (rtot (alignedRk M r)) =>
            ABlockW (optWG β) β (alignedR M r) t.1 t.2) := fun t => hgram t.1 t.2
      have h := TendstoInProbPi.comp_continuous
        (continuous_eigVal_symMat (p := rtot (alignedRk M r)) (j : ℕ)).continuousAt hpi
      rw [eigVal_symMat (isHermitian_ABlockW (optWG β) β (alignedR M r)), hev, hSz,
        add_zero] at h
      refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
      exact eigVal_symMat (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ)
    have hy : TendstoInProb μ (fun N ω => ‖m.VtVWGcol (optWG β) N ω j‖ ^ 2) 0 := by
      have hpi : TendstoInProbPi μ
          (fun N ω => fun q : Fin (rtot (alignedRk M r)) => m.VtVWG (optWG β) N ω q j)
          (fun q : Fin (rtot (alignedRk M r)) =>
            BBlockW (optWG β) β (alignedR M r) q j) := fun q => hVtV q j
      have hcont : Continuous fun z : Fin (rtot (alignedRk M r)) → ℝ => ∑ q, z q ^ 2 :=
        continuous_finsetSum _ fun q _ => (continuous_apply q).pow 2
      have h := TendstoInProbPi.comp_continuous hcont.continuousAt hpi
      have hlim : ∑ q : Fin (rtot (alignedRk M r)),
          BBlockW (optWG β) β (alignedR M r) q j ^ 2 = 0 := by
        have hc0 : ‖BBlockWcol (optWG β) β (alignedR M r) j‖ ^ 2 = 0 := by
          rw [norm_sq_BBlockWcol_aligned β hsq j, hSz]
        rw [← hc0]
        exact (normSq_toLp fun q => BBlockW (optWG β) β (alignedR M r) q j).symm
      rw [hlim] at h
      refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
      exact (normSq_toLp fun q => m.VtVWG (optWG β) N ω q j).symm
    have hnum : TendstoInProb μ
        (fun N ω => ‖specProjIdx (m.gramWG (optWG β) N ω)
          (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ)
          (m.VtVWGcol (optWG β) N ω j)‖ ^ 2) 0 := by
      refine TendstoInProb.of_le (fun N => Filter.Eventually.of_forall fun ω => ?_) hy
      have hle : ‖specProjIdx (m.gramWG (optWG β) N ω)
          (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ)
          (m.VtVWGcol (optWG β) N ω j)‖ ≤ ‖m.VtVWGcol (optWG β) N ω j‖ :=
        norm_specProj_le _ _ _
      rw [sub_zero, abs_of_nonneg (sq_nonneg _)]
      nlinarith [hle, norm_nonneg (specProjIdx (m.gramWG (optWG β) N ω)
        (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ) (m.VtVWGcol (optWG β) N ω j))]
    have hdiv := TendstoInProb.div hnum hden one_ne_zero
    rw [zero_div] at hdiv
    exact hdiv
  · -- branch P: `0 < S_j`, the two gaps hold and the quotient functional is continuous
    have hgapj := topGap_ABlockW_optWG_aligned_of_sep β h0 h1 hsep hM j hpos
    have hgapj1 := topGap_succ_ABlockW_optWG_aligned_of_sep β h0 h1 hsep hM j hpos
    have hposEV : 0 < eigVal (ABlockW (optWG β) β (alignedR M r))
        (isHermitian_ABlockW (optWG β) β (alignedR M r)) (j : ℕ) := by
      rw [hev]; linarith
    have hconv : TendstoInProbPi μ
        (fun N ω => Sum.elim
          (fun t : Fin (rtot (alignedRk M r)) × Fin (rtot (alignedRk M r)) =>
            m.gramWG (optWG β) N ω t.1 t.2)
          (fun t : Fin (rtot (alignedRk M r)) × Fin r => m.VtVWG (optWG β) N ω t.1 t.2))
        (Sum.elim
          (fun t : Fin (rtot (alignedRk M r)) × Fin (rtot (alignedRk M r)) =>
            ABlockW (optWG β) β (alignedR M r) t.1 t.2)
          (fun t : Fin (rtot (alignedRk M r)) × Fin r =>
            BBlockW (optWG β) β (alignedR M r) t.1 t.2)) := by
      rintro (t | t)
      · exact hgram t.1 t.2
      · exact hVtV t.1 t.2
    have hall := TendstoInProbPi.comp_continuous
      (continuousAt_compFun (isHermitian_ABlockW (optWG β) β (alignedR M r)) hgapj hgapj1
        hposEV (BBlockW (optWG β) β (alignedR M r)) j) hconv
    have hval : ‖specProjIdx (ABlockW (optWG β) β (alignedR M r))
          (isHermitian_ABlockW (optWG β) β (alignedR M r)) (j : ℕ)
          (WithLp.toLp 2 fun q => BBlockW (optWG β) β (alignedR M r) q j)‖ ^ 2
        / eigVal (ABlockW (optWG β) β (alignedR M r))
          (isHermitian_ABlockW (optWG β) β (alignedR M r)) (j : ℕ)
        = Sagg β j / (Sagg β j + 1) := by
      rw [hev]
      exact compLimit_aligned_of_sep β h0 h1 hsep hM j hpos
    rw [compFun_eq (isHermitian_ABlockW (optWG β) β (alignedR M r))
      (BBlockW (optWG β) β (alignedR M r)) (j : ℕ) j, hval] at hall
    refine hall.congr fun N => Filter.Eventually.of_forall fun ω => ?_
    exact compFun_eq (m.isHermitian_gramWG (optWG β) N ω) (m.VtVWG (optWG β) N ω) (j : ℕ) j

theorem thm_rank_r_svdstack_component (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hS : StrictAnti (Sagg β))
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) (j : Fin r) :
    TendstoInProb μ (fun N ω => m.compRGW (optWG β) j N ω) (Sagg β j / (Sagg β j + 1)) :=
  m.thm_rank_r_svdstack_component_of_sep c β hc hβdef hR hM (Sagg_sep_of_strictAnti hS) law hG j

/-- **`thm:rank_r_svdstack`, component clause at the canonical frame** (`main_paper.tex:2405`,
Track D item D4). `v̂_j` is column `j` of `V̂_svdstack(W⋆) = Ṽ_{W⋆}ᵀ Q_r Λ_r^{-1/2}` at the
canonical top-`r` eigenframe `(Q_r, Λ_r) = (topEigMat, topEigVal)` of `Ṽ_{W⋆} Ṽ_{W⋆}ᵀ`, whose
columns are ordered by decreasing eigenvalue by construction. No frame hypothesis appears.

Route: the frame entry squared is `⟪q_j, y_j⟫² / λ_j` (`vhatSvdstackGW_entry_sq`), and at the
canonical frame `q_j = vEig` and `λ_j = eigVal`. At `0 < S_j` the sorted index `j` of the
random Gram matrix is simple with probability tending to one, and there the quantity is
`compRGW`; at `S_j = 0` the numerator is squeezed by `‖y_j‖² → 0`. -/
theorem thm_rank_r_svdstack_component_eig_of_sep
    (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hsep : SaggSep β)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) (j : Fin r) :
    TendstoInProb μ
      (fun N ω => (((m.vhatSvdstackGW (optWG β) N ω
          (topEigMat (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM))
          (topEigVal (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM)))ᵀ * m.V N) j j) ^ 2)
      (Sagg β j / (Sagg β j + 1)) := by
  have hβ01 : ∀ i k, 0 ≤ β i k ∧ β i k < 1 := fun i k => by
    rw [hβdef i k]
    exact beta_mem_Ico (hc i)
  have h0 : ∀ i k, 0 ≤ β i k := fun i k => (hβ01 i k).1
  have h1 : ∀ i k, β i k < 1 := fun i k => (hβ01 i k).2
  have hsq : ∀ i k, β i k ^ 2 < 1 := fun i k => by nlinarith [h0 i k, h1 i k]
  have hrp : r ≤ Fintype.card (Fin (rtot (alignedRk M r))) := by
    simpa using le_rtot_alignedRk (M := M) (r := r) hM
  have hjc : (j : ℕ) < Fintype.card (Fin (rtot (alignedRk M r))) :=
    lt_card_rtot_alignedRk hM j.isLt
  have hev := eigVal_ABlockW_aligned β h0 h1 hsep.1 hM j
  have hgram := m.gramWG_tendsto_aligned c β hβdef hR law hG
  have hVtV := m.VtVWG_tendsto_aligned c β hβdef hR law
  have hLnn : ∀ (N : ℕ) (ω : Ω N),
      0 ≤ eigVal (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ) := by
    intro N ω
    rw [← topEigVal_eq_eigVal (m.isHermitian_gramWG (optWG β) N ω) hrp j]
    exact topEigVal_nonneg hrp (m.posSemidef_gramWG (optWG β) N ω) j
  have hentry : ∀ (N : ℕ) (ω : Ω N),
      (((m.vhatSvdstackGW (optWG β) N ω (topEigMat (m.isHermitian_gramWG (optWG β) N ω) hrp)
            (topEigVal (m.isHermitian_gramWG (optWG β) N ω) hrp))ᵀ * m.V N) j j) ^ 2
        = ⟪vEig (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ),
            m.VtVWGcol (optWG β) N ω j⟫_ℝ ^ 2
          / eigVal (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ) := by
    intro N ω
    rw [m.vhatSvdstackGW_entry_sq (optWG β) N ω _ _ j
        (topEigVal_nonneg hrp (m.posSemidef_gramWG (optWG β) N ω) j),
      frameCol_topEigMat_eq_vEig, topEigVal_eq_eigVal]
  suffices h : TendstoInProb μ
      (fun N ω => ⟪vEig (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ),
          m.VtVWGcol (optWG β) N ω j⟫_ℝ ^ 2
        / eigVal (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ))
      (Sagg β j / (Sagg β j + 1)) from
    h.congr fun N => Filter.Eventually.of_forall fun ω => (hentry N ω).symm
  rcases eq_or_lt_of_le (Sagg_nonneg hsq j) with hzero | hpos
  · -- branch Z: no gap, no simplicity; Cauchy-Schwarz and the squeeze
    have hSz : Sagg β j = 0 := hzero.symm
    have hlimval : Sagg β j / (Sagg β j + 1) = 0 := by rw [hSz]; norm_num
    rw [hlimval]
    have hden : TendstoInProb μ
        (fun N ω => eigVal (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω)
          (j : ℕ)) 1 := by
      have hpi : TendstoInProbPi μ
          (fun N ω => fun t : Fin (rtot (alignedRk M r)) × Fin (rtot (alignedRk M r)) =>
            m.gramWG (optWG β) N ω t.1 t.2)
          (fun t : Fin (rtot (alignedRk M r)) × Fin (rtot (alignedRk M r)) =>
            ABlockW (optWG β) β (alignedR M r) t.1 t.2) := fun t => hgram t.1 t.2
      have h := TendstoInProbPi.comp_continuous
        (continuous_eigVal_symMat (p := rtot (alignedRk M r)) (j : ℕ)).continuousAt hpi
      rw [eigVal_symMat (isHermitian_ABlockW (optWG β) β (alignedR M r)), hev, hSz,
        add_zero] at h
      refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
      exact eigVal_symMat (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ)
    have hy : TendstoInProb μ (fun N ω => ‖m.VtVWGcol (optWG β) N ω j‖ ^ 2) 0 := by
      have hpi : TendstoInProbPi μ
          (fun N ω => fun q : Fin (rtot (alignedRk M r)) => m.VtVWG (optWG β) N ω q j)
          (fun q : Fin (rtot (alignedRk M r)) =>
            BBlockW (optWG β) β (alignedR M r) q j) := fun q => hVtV q j
      have hcont : Continuous fun z : Fin (rtot (alignedRk M r)) → ℝ => ∑ q, z q ^ 2 :=
        continuous_finsetSum _ fun q _ => (continuous_apply q).pow 2
      have h := TendstoInProbPi.comp_continuous hcont.continuousAt hpi
      have hlim : ∑ q : Fin (rtot (alignedRk M r)),
          BBlockW (optWG β) β (alignedR M r) q j ^ 2 = 0 := by
        have hc0 : ‖BBlockWcol (optWG β) β (alignedR M r) j‖ ^ 2 = 0 := by
          rw [norm_sq_BBlockWcol_aligned β hsq j, hSz]
        rw [← hc0]
        exact (normSq_toLp fun q => BBlockW (optWG β) β (alignedR M r) q j).symm
      rw [hlim] at h
      refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
      exact (normSq_toLp fun q => m.VtVWG (optWG β) N ω q j).symm
    have hZ := TendstoInProb.div hy hden one_ne_zero
    rw [zero_div] at hZ
    refine TendstoInProb.of_le (fun N => Filter.Eventually.of_forall fun ω => ?_) hZ
    rw [sub_zero]
    exact inner_sq_div_le _ _ (norm_vEig _ _ hjc) (hLnn N ω)
  · -- branch P: the sorted index `j` is simple with probability tending to one
    have hgapj := topGap_ABlockW_optWG_aligned_of_sep β h0 h1 hsep hM j hpos
    have hgapj1 := topGap_succ_ABlockW_optWG_aligned_of_sep β h0 h1 hsep hM j hpos
    have hwhp := specIdxSimple_whp_of_tendsto
      (isHermitian_ABlockW (optWG β) β (alignedR M r))
      (fun N ω => m.isHermitian_gramWG (optWG β) N ω) (fun p q => hgram p q) hjc hgapj hgapj1
    refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_
      (m.thm_rank_r_svdstack_component_of_sep c β hc hβdef hR hM hsep law hG j)
    refine tendsto_measure_zero_of_subset (t := fun N =>
      {ω | ¬ SimpleIdx (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ)})
      (fun N ω hω => ?_) hwhp
    by_contra hnot
    simp only [Set.mem_ofPred_eq, not_not] at hnot
    refine hω ?_
    show _ = m.compRGW (optWG β) j N ω
    rw [UnalignedModelR.compRGW, normSq_specProjIdx_eq_inner_sq hnot]

theorem thm_rank_r_svdstack_component_eig (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hS : StrictAnti (Sagg β))
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) (j : Fin r) :
    TendstoInProb μ
      (fun N ω => (((m.vhatSvdstackGW (optWG β) N ω
          (topEigMat (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM))
          (topEigVal (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM)))ᵀ * m.V N) j j) ^ 2)
      (Sagg β j / (Sagg β j + 1)) :=
  m.thm_rank_r_svdstack_component_eig_of_sep c β hc hβdef hR hM
    (Sagg_sep_of_strictAnti hS) law hG j

theorem thm_rank_r_svdstack_component_of_model
    (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise) (j : Fin r) :
    TendstoInProb μ
      (fun N ω => (((m.vhatSvdstackGW (optWG β) N ω
          (topEigMat (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM))
          (topEigVal (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM)))ᵀ * m.V N) j j) ^ 2)
      (Sagg β j / (Sagg β j + 1)) :=
  m.thm_rank_r_svdstack_component_eig_of_sep c β hc hβdef hR hM
    (m.saggSep_of_model c β hc hβdef) law hG j

/-- **`thm:rank_r_svdstack`, component clause at any measurable sorted frame**
(`main_paper.tex:2405`, Track D item D4). `IsTopEigFrame` alone does not order the columns of
`Q`, so a component-labelled statement needs the extra conjunct `lam N ω k = λ_k` in the good
event (modeling choice M3 of `notes/archive/rankr_D4_plan.md`). The canonical frame satisfies both
conjuncts with no hypothesis; that is
`thm_rank_r_svdstack_component_eig_of_sep`. -/
theorem thm_rank_r_svdstack_component_frame_of_sep
    (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hsep : SaggSep β)
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise)
    (Q : ∀ N, Ω N → Matrix (Fin (rtot (alignedRk M r))) (Fin r) ℝ)
    (lam : ∀ N, Ω N → Fin r → ℝ)
    (hQ : Tendsto (fun N => μ N {ω | ¬ (IsTopEigFrame (m.gramWG (optWG β) N ω)
      (m.isHermitian_gramWG (optWG β) N ω) (Q N ω) (lam N ω) ∧
      ∀ k, lam N ω k = eigVal (m.gramWG (optWG β) N ω)
        (m.isHermitian_gramWG (optWG β) N ω) (k : ℕ))}) atTop (𝓝 0)) (j : Fin r) :
    TendstoInProb μ
      (fun N ω => (((m.vhatSvdstackGW (optWG β) N ω (Q N ω) (lam N ω))ᵀ * m.V N) j j) ^ 2)
      (Sagg β j / (Sagg β j + 1)) := by
  classical
  have hβ01 : ∀ i k, 0 ≤ β i k ∧ β i k < 1 := fun i k => by
    rw [hβdef i k]
    exact beta_mem_Ico (hc i)
  have h0 : ∀ i k, 0 ≤ β i k := fun i k => (hβ01 i k).1
  have h1 : ∀ i k, β i k < 1 := fun i k => (hβ01 i k).2
  have hsq : ∀ i k, β i k ^ 2 < 1 := fun i k => by nlinarith [h0 i k, h1 i k]
  have hrp : r ≤ Fintype.card (Fin (rtot (alignedRk M r))) := by
    simpa using le_rtot_alignedRk (M := M) (r := r) hM
  have hjc : (j : ℕ) < Fintype.card (Fin (rtot (alignedRk M r))) :=
    lt_card_rtot_alignedRk hM j.isLt
  have hev := eigVal_ABlockW_aligned β h0 h1 hsep.1 hM j
  have hgram := m.gramWG_tendsto_aligned c β hβdef hR law hG
  have hVtV := m.VtVWG_tendsto_aligned c β hβdef hR law
  have hLnn : ∀ (N : ℕ) (ω : Ω N),
      0 ≤ eigVal (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ) := by
    intro N ω
    rw [← topEigVal_eq_eigVal (m.isHermitian_gramWG (optWG β) N ω) hrp j]
    exact topEigVal_nonneg hrp (m.posSemidef_gramWG (optWG β) N ω) j
  -- the entry identity, on the frame event
  have hframe : ∀ (N : ℕ) (ω : Ω N),
      IsTopEigFrame (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω)
        (Q N ω) (lam N ω) →
      (∀ k, lam N ω k = eigVal (m.gramWG (optWG β) N ω)
        (m.isHermitian_gramWG (optWG β) N ω) (k : ℕ)) →
      (((m.vhatSvdstackGW (optWG β) N ω (Q N ω) (lam N ω))ᵀ * m.V N) j j) ^ 2
        = ⟪frameCol (Q N ω) j, m.VtVWGcol (optWG β) N ω j⟫_ℝ ^ 2
          / eigVal (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ) := by
    intro N ω _ hQ2
    rw [m.vhatSvdstackGW_entry_sq (optWG β) N ω _ _ j (by rw [hQ2 j]; exact hLnn N ω), hQ2 j]
  rcases eq_or_lt_of_le (Sagg_nonneg hsq j) with hzero | hpos
  · -- branch Z: the frame column is a unit vector on the event, so Cauchy-Schwarz squeezes
    have hSz : Sagg β j = 0 := hzero.symm
    have hlimval : Sagg β j / (Sagg β j + 1) = 0 := by rw [hSz]; norm_num
    rw [hlimval]
    have hden : TendstoInProb μ
        (fun N ω => eigVal (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω)
          (j : ℕ)) 1 := by
      have hpi : TendstoInProbPi μ
          (fun N ω => fun t : Fin (rtot (alignedRk M r)) × Fin (rtot (alignedRk M r)) =>
            m.gramWG (optWG β) N ω t.1 t.2)
          (fun t : Fin (rtot (alignedRk M r)) × Fin (rtot (alignedRk M r)) =>
            ABlockW (optWG β) β (alignedR M r) t.1 t.2) := fun t => hgram t.1 t.2
      have h := TendstoInProbPi.comp_continuous
        (continuous_eigVal_symMat (p := rtot (alignedRk M r)) (j : ℕ)).continuousAt hpi
      rw [eigVal_symMat (isHermitian_ABlockW (optWG β) β (alignedR M r)), hev, hSz,
        add_zero] at h
      refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
      exact eigVal_symMat (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ)
    have hy : TendstoInProb μ (fun N ω => ‖m.VtVWGcol (optWG β) N ω j‖ ^ 2) 0 := by
      have hpi : TendstoInProbPi μ
          (fun N ω => fun q : Fin (rtot (alignedRk M r)) => m.VtVWG (optWG β) N ω q j)
          (fun q : Fin (rtot (alignedRk M r)) =>
            BBlockW (optWG β) β (alignedR M r) q j) := fun q => hVtV q j
      have hcont : Continuous fun z : Fin (rtot (alignedRk M r)) → ℝ => ∑ q, z q ^ 2 :=
        continuous_finsetSum _ fun q _ => (continuous_apply q).pow 2
      have h := TendstoInProbPi.comp_continuous hcont.continuousAt hpi
      have hlim : ∑ q : Fin (rtot (alignedRk M r)),
          BBlockW (optWG β) β (alignedR M r) q j ^ 2 = 0 := by
        have hc0 : ‖BBlockWcol (optWG β) β (alignedR M r) j‖ ^ 2 = 0 := by
          rw [norm_sq_BBlockWcol_aligned β hsq j, hSz]
        rw [← hc0]
        exact (normSq_toLp fun q => BBlockW (optWG β) β (alignedR M r) q j).symm
      rw [hlim] at h
      refine h.congr fun N => Filter.Eventually.of_forall fun ω => ?_
      exact (normSq_toLp fun q => m.VtVWG (optWG β) N ω q j).symm
    have hZ := TendstoInProb.div hy hden one_ne_zero
    rw [zero_div] at hZ
    have hg : TendstoInProb μ
        (fun N ω => if (IsTopEigFrame (m.gramWG (optWG β) N ω)
            (m.isHermitian_gramWG (optWG β) N ω) (Q N ω) (lam N ω) ∧
            ∀ k, lam N ω k = eigVal (m.gramWG (optWG β) N ω)
              (m.isHermitian_gramWG (optWG β) N ω) (k : ℕ))
          then (((m.vhatSvdstackGW (optWG β) N ω (Q N ω) (lam N ω))ᵀ * m.V N) j j) ^ 2
          else 0) 0 := by
      refine TendstoInProb.of_le (fun N => Filter.Eventually.of_forall fun ω => ?_) hZ
      rw [sub_zero]
      by_cases hgood : (IsTopEigFrame (m.gramWG (optWG β) N ω)
          (m.isHermitian_gramWG (optWG β) N ω) (Q N ω) (lam N ω) ∧
          ∀ k, lam N ω k = eigVal (m.gramWG (optWG β) N ω)
            (m.isHermitian_gramWG (optWG β) N ω) (k : ℕ))
      · rw [if_pos hgood, hframe N ω hgood.1 hgood.2]
        exact inner_sq_div_le _ _ (norm_frameCol_of_ortho hgood.1.frame.ortho j) (hLnn N ω)
      · rw [if_neg hgood, abs_zero]
        exact div_nonneg (sq_nonneg _) (hLnn N ω)
    refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_ hg
    refine tendsto_measure_zero_of_subset (t := fun N => {ω | ¬ (IsTopEigFrame
      (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (Q N ω) (lam N ω) ∧
      ∀ k, lam N ω k = eigVal (m.gramWG (optWG β) N ω)
        (m.isHermitian_gramWG (optWG β) N ω) (k : ℕ))}) (fun N ω hω => ?_) hQ
    by_contra hnot
    simp only [Set.mem_ofPred_eq, not_not] at hnot
    exact hω (if_pos hnot).symm
  · -- branch P: the frame column is the sorted eigenvector, up to a sign, where `j` is simple
    have hgapj := topGap_ABlockW_optWG_aligned_of_sep β h0 h1 hsep hM j hpos
    have hgapj1 := topGap_succ_ABlockW_optWG_aligned_of_sep β h0 h1 hsep hM j hpos
    have hwhp := specIdxSimple_whp_of_tendsto
      (isHermitian_ABlockW (optWG β) β (alignedR M r))
      (fun N ω => m.isHermitian_gramWG (optWG β) N ω) (fun p q => hgram p q) hjc hgapj hgapj1
    refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_
      (m.thm_rank_r_svdstack_component_of_sep c β hc hβdef hR hM hsep law hG j)
    refine tendsto_measure_zero_of_subset (t := fun N =>
      {ω | ¬ (IsTopEigFrame (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω)
        (Q N ω) (lam N ω) ∧ ∀ k, lam N ω k = eigVal (m.gramWG (optWG β) N ω)
          (m.isHermitian_gramWG (optWG β) N ω) (k : ℕ))} ∪
      {ω | ¬ SimpleIdx (m.gramWG (optWG β) N ω) (m.isHermitian_gramWG (optWG β) N ω) (j : ℕ)})
      (fun N ω hω => ?_) (tendsto_measure_zero_union hQ hwhp)
    by_contra hnot
    simp only [Set.mem_union, Set.mem_ofPred_eq, not_or, not_not] at hnot
    obtain ⟨hgood, hsimpleω⟩ := hnot
    have hlamj : lam N ω j
        = (m.isHermitian_gramWG (optWG β) N ω).eigenvalues₀ ⟨(j : ℕ), hjc⟩ := by
      rw [hgood.2 j, eigVal_eq]
    have hq : toOp (m.gramWG (optWG β) N ω) (frameCol (Q N ω) j)
        = (m.isHermitian_gramWG (optWG β) N ω).eigenvalues₀ ⟨(j : ℕ), hjc⟩ •
          frameCol (Q N ω) j := by
      rw [← hlamj]
      exact toOp_frameCol hgood.1.eig j
    refine hω ?_
    show _ = m.compRGW (optWG β) j N ω
    rw [hframe N ω hgood.1 hgood.2, UnalignedModelR.compRGW,
      inner_sq_eq_normSq_specProjIdx hsimpleω hjc hq
        (norm_frameCol_of_ortho hgood.1.frame.ortho j)]

theorem thm_rank_r_svdstack_component_frame (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hS : StrictAnti (Sagg β))
    (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise)
    (Q : ∀ N, Ω N → Matrix (Fin (rtot (alignedRk M r))) (Fin r) ℝ)
    (lam : ∀ N, Ω N → Fin r → ℝ)
    (hQ : Tendsto (fun N => μ N {ω | ¬ (IsTopEigFrame (m.gramWG (optWG β) N ω)
      (m.isHermitian_gramWG (optWG β) N ω) (Q N ω) (lam N ω) ∧
      ∀ k, lam N ω k = eigVal (m.gramWG (optWG β) N ω)
        (m.isHermitian_gramWG (optWG β) N ω) (k : ℕ))}) atTop (𝓝 0)) (j : Fin r) :
    TendstoInProb μ
      (fun N ω => (((m.vhatSvdstackGW (optWG β) N ω (Q N ω) (lam N ω))ᵀ * m.V N) j j) ^ 2)
      (Sagg β j / (Sagg β j + 1)) :=
  m.thm_rank_r_svdstack_component_frame_of_sep c β hc hβdef hR hM
    (Sagg_sep_of_strictAnti hS) law hG Q lam hQ j

end UnalignedModelR

end StackedSVD
