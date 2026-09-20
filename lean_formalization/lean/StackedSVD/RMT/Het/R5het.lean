/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.R4het
import StackedSVD.RMT.Het.MPhet

/-!
# Item H4: the heteroscedastic analytic core

Task H4 of `notes/archive/plan_heterolaw_A.md` (sections 2.4, 3.1, 3.2, 3.7 and the H4 row of
section 4). Nothing here is Gaussian and nothing here uses independence: the file turns the
interface `ResolventLimitsHet` of `RMT/Het/R4het.lean` (an edge bound and six resolvent
limits at every real `z > b`) into the three conclusions

* `lamMax_tendstoInProb_het`: `λ_max(X_wᵀ X_w) → ρ` in probability;
* `tendsto_measure_topSimple_het`: the top eigenvalue of `X_wᵀ X_w` is simple with
  probability tending to `1`;
* `align_tendstoInProb_het`: `overlap(X_w, v) → L(w)` in probability.

The route is `RMT/R5.lean` with the column split of `RMT/Het/Split.lean` in place of the
row split. `X_w X_wᵀ = W₀' + q qᵀ` (`gram_eq_het`) is a rank-one positive update on the
`n` side, so `R4.lamMax_eq`, `R4.topSimple` and `R4.qform2_antitoneOn` apply verbatim
there; `Het.lamMax_gram_comm` and the new `Het.topSimple_transpose_mul_self` move the
eigenvalue and its simplicity to the `d` side, and `stackPerfW_eq_inv_lam_qform2`
(`RMT/Het/Duality.lean`) replaces R5's numerator bracket by the single identity
`overlap = 1/(λ · qᵀ G₀'(λ)² q)`.

## Choices

1. The three conclusions stay model-free, over the abstract scalar facts `HetScalarFacts`
   of `R4het.lean`, exactly as R5 is stated over `ResolventLimits`. The concrete branch
   enters through `hetScalarFacts_of_assumption4`: under `Assumption4`, `0 < c`, and
   `bHet ≤ b < rhoHet`, the closed forms `Phihet`, `Psihet`, `PhihetDeriv`, `PsihetDeriv`
   of `RMT/Het/MPhet.lean` satisfy every field. The corollaries `_of_assumption4` are the
   statements of the H4 row of the plan.
2. `HetScalarFacts` has no sign of `b`; `b_nonneg_of_edge` recovers `0 ≤ b` from the edge
   field, because `W₀'` is positive semidefinite. This is what makes `ρ > 0` and the
   denominator bracket of the overlap positive.
3. The three statements carry `[∀ N, IsProbabilityMeasure (μ N)]` (audit item F1,
   `notes/archive/audit_het_skeleton_2026-08-31.md`); `R4het.lean` lacked it.
4. The measurability of the edge event (audit item F2) is proved here, not in
   `Split.lean`: `W₀' = Y Yᵀ` with `Y = Σ^{1/2} E⊥`, so `λ_max(W₀') = gramLamMax Yᵀ`.

Numeric check (a session script, `check_h4.py`, not kept; seed 20260901, `d = 400`,
`M = 3`):
`1 + Fhet(ρ) = 1.1e-16`; the secular bracket at `ρ ∓ 0.5` is `-0.018 / +0.076`;
`overlap = 1/(λ · qform2) = 0.63734` both ways, `Lw = 0.63927`; `λ_max = 15.40` vs
`ρ = 15.73`.

Paper: `main_paper.tex` lines 1409 to 1440 (`thm:stacksvd_weighted`, `eq:assumption4`).
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

/-! ### Simplicity moves between the two Gram matrices -/

namespace Het

variable {p q : ℕ}

/-- The top eigenvalue of a positive semidefinite matrix is nonnegative. -/
theorem lamMax_nonneg_of_posSemidef {A : Matrix (Fin p) (Fin p) ℝ}
    (hA : A.IsHermitian) (hpsd : A.PosSemidef) (hp : 0 < p) : 0 ≤ lamMax A hA := by
  obtain ⟨j, hj⟩ := R4.exists_eigenvalues_eq_lamMax hA hp
  rw [← hj]
  exact hpsd.eigenvalues_nonneg j

set_option linter.deprecated false in
/-- `v ↦ X v` embeds the `t`-eigenspace of `Xᵀ X` into the `t`-eigenspace of `X Xᵀ` when
`t ≠ 0`, so the first has the smaller dimension. -/
theorem finrank_eigenspace_gram_le (X : Matrix (Fin p) (Fin q) ℝ) {t : ℝ} (ht : t ≠ 0) :
    Module.finrank ℝ (Module.End.eigenspace (toOp (Xᵀ * X)) t)
      ≤ Module.finrank ℝ (Module.End.eigenspace (toOp (X * Xᵀ)) t) := by
  set Ed := Module.End.eigenspace (toOp (Xᵀ * X)) t with hEd
  set En := Module.End.eigenspace (toOp (X * Xᵀ)) t with hEn
  have hmem : ∀ v ∈ Ed, Matrix.toEuclideanLin X v ∈ En := by
    intro v hv
    rw [hEd, R4.mem_eigenspace_iff'] at hv
    rw [hEn, R4.mem_eigenspace_iff']
    simp only [Matrix.toEuclideanLin_apply, WithLp.ofLp_toLp]
    rw [← Matrix.mulVec_mulVec, Matrix.mulVec_mulVec (WithLp.ofLp v) Xᵀ X, hv,
      Matrix.mulVec_smul]
  let f : Ed →ₗ[ℝ] En := (Matrix.toEuclideanLin X).restrict hmem
  have hinj : Function.Injective f := by
    rw [← LinearMap.ker_eq_bot, LinearMap.ker_eq_bot']
    intro v hv
    have h1 : Matrix.toEuclideanLin X v = 0 := by
      have := congrArg Subtype.val hv
      simpa [f, LinearMap.restrict_apply] using this
    have h2 : X *ᵥ WithLp.ofLp (v : EuclideanSpace ℝ (Fin q)) = 0 := by
      have := congrArg WithLp.ofLp h1
      simpa [Matrix.toEuclideanLin_apply] using this
    have hv' := (R4.mem_eigenspace_iff' _ _ _).mp v.2
    rw [← Matrix.mulVec_mulVec, h2, Matrix.mulVec_zero] at hv'
    have h3 : WithLp.ofLp (v : EuclideanSpace ℝ (Fin q)) = 0 :=
      (smul_eq_zero.mp hv'.symm).resolve_left ht
    apply Subtype.ext
    exact WithLp.ofLp_injective 2 (by simpa using h3)
  exact LinearMap.finrank_le_finrank_of_injective hinj

/-- **Simplicity moves from `X Xᵀ` to `Xᵀ X`** when the shared top eigenvalue is positive:
`v ↦ X v` embeds the top eigenspace of `Xᵀ X` into that of `X Xᵀ`, and the former is not
`⊥`. Step 4 of `Het.topProj_transpose_eq` in the direction that item H4 needs. -/
theorem topSimple_transpose_mul_self (hp : 0 < p) (hq : 0 < q) (X : Matrix (Fin p) (Fin q) ℝ)
    (hlam : 0 < gramLamMax X)
    (hsimple : TopSimple (X * Xᵀ) (isHermitian_mul_transpose_self X)) :
    TopSimple (Xᵀ * X) (isHermitian_transpose_mul_self X) := by
  have hA := isHermitian_transpose_mul_self X
  rw [TopSimple, topSpace_eq_eigenspace] at hsimple ⊢
  rw [← lamMax_gram_comm hp hq X] at hsimple
  have hle := finrank_eigenspace_gram_le X hlam.ne'
  have hpos : 0 < Module.finrank ℝ (Module.End.eigenspace (toOp (Xᵀ * X)) (gramLamMax X)) := by
    obtain ⟨j, hj⟩ := R4.exists_eigenvalues_eq_lamMax hA hq
    have hmem : hA.eigenvectorBasis j
        ∈ Module.End.eigenspace (toOp (Xᵀ * X)) (gramLamMax X) := by
      rw [R4.mem_eigenspace_iff']
      change (Xᵀ * X) *ᵥ ⇑(hA.eigenvectorBasis j) = gramLamMax X • ⇑(hA.eigenvectorBasis j)
      rw [show gramLamMax X = hA.eigenvalues j from hj.symm]
      exact hA.mulVec_eigenvectorBasis j
    have hne : hA.eigenvectorBasis j ≠ 0 := hA.eigenvectorBasis.orthonormal.ne_zero j
    refine Nat.pos_of_ne_zero fun h0 => hne ?_
    have hbot := Submodule.finrank_eq_zero.mp h0
    rw [hbot] at hmem
    exact (Submodule.mem_bot ℝ).mp hmem
  change Module.finrank ℝ (Module.End.eigenspace (toOp (Xᵀ * X)) (gramLamMax X)) = 1
  omega

end Het

/-! ### The concrete scalar facts (audit item F3)

`HetScalarFacts` of `R4het.lean` lists what the analytic core needs from the limiting
secular function. `RMT/Het/MPhet.lean` proves the two identities (`one_add_F_eq_zero_iff`,
`overlap_identity_het`); the sign change, the continuity and the positivity of the
denominator are proved here from `one_add_F_eq`, `sPhys_strictMonoOn` and
`rhoHet_mul_FhetDeriv`. -/

namespace MPhet

open Scalars

variable {M : ℕ} {θ c w : Fin M → ℝ}

theorem continuousOn_ghet (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) (i : Fin M) :
    ContinuousOn (ghet c w i) (Set.Ioi (bHet c w)) := by
  intro z hz
  have hz0 : z ≠ 0 := (lt_trans (bHet_pos hc hw) hz).ne'
  have hden : z * (1 + w i ^ 2 * sPhys c w z) ≠ 0 :=
    mul_ne_zero hz0 (one_add_mul_sPhys_pos hc hw hz i).ne'
  have hs := continuousAt_sPhys hc hw hz
  refine ContinuousAt.continuousWithinAt ?_
  change ContinuousAt (fun z => -1 / (z * (1 + w i ^ 2 * sPhys c w z))) z
  exact continuousAt_const.div
    (continuousAt_id.mul (continuousAt_const.add (continuousAt_const.mul hs))) hden

theorem continuousOn_Phihet (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    ContinuousOn (Phihet θ c w) (Set.Ioi (bHet c w)) := by
  unfold Phihet
  exact continuousOn_finsetSum _ fun i _ => continuousOn_const.mul (continuousOn_ghet hc hw i)

theorem continuousOn_Psihet (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    ContinuousOn (Psihet c w) (Set.Ioi (bHet c w)) := by
  unfold Psihet
  exact continuousOn_finsetSum _ fun i _ => continuousOn_const.mul (continuousOn_ghet hc hw i)

theorem continuousOn_sPhysDeriv (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    ContinuousOn (sPhysDeriv c w) (Set.Ioi (bHet c w)) := by
  intro z hz
  refine ContinuousAt.continuousWithinAt ?_
  have hs := continuousAt_sPhys hc hw hz
  have hmem : sPhys c w z ∈ Set.Ioo (sLo w) 0 := ⟨sLo_lt_sPhys hc hw hz, sPhys_neg hc hw hz⟩
  have hz' : ContinuousAt (zfunDeriv c w) (sPhys c w z) :=
    (continuousOn_zfunDeriv hw).continuousAt (Ioo_mem_nhds hmem.1 hmem.2)
  have hne : zfunDeriv c w (sPhys c w z) ≠ 0 :=
    (zfunDeriv_pos hc hw (sStar_lt_sPhys hc hw hz) (sPhys_neg hc hw hz)).ne'
  change ContinuousAt (fun z => (zfunDeriv c w (sPhys c w z))⁻¹) z
  exact (hz'.comp hs).inv₀ hne

theorem continuousOn_ghetDeriv (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) (i : Fin M) :
    ContinuousOn (ghetDeriv c w i) (Set.Ioi (bHet c w)) := by
  intro z hz
  have hz0 : z ≠ 0 := (lt_trans (bHet_pos hc hw) hz).ne'
  have hden : (z * (1 + w i ^ 2 * sPhys c w z)) ^ 2 ≠ 0 :=
    pow_ne_zero 2 (mul_ne_zero hz0 (one_add_mul_sPhys_pos hc hw hz i).ne')
  have hs := continuousAt_sPhys hc hw hz
  have hs' := (continuousOn_sPhysDeriv hc hw).continuousAt (Ioi_mem_nhds hz)
  refine ContinuousAt.continuousWithinAt ?_
  change ContinuousAt (fun z => (1 + w i ^ 2 * sPhys c w z + z * w i ^ 2 * sPhysDeriv c w z)
    / (z * (1 + w i ^ 2 * sPhys c w z)) ^ 2) z
  refine ContinuousAt.div ?_ ?_ hden
  · exact (continuousAt_const.add (continuousAt_const.mul hs)).add
      ((continuousAt_id.mul continuousAt_const).mul hs')
  · exact (continuousAt_id.mul (continuousAt_const.add (continuousAt_const.mul hs))).pow 2

theorem continuousOn_PhihetDeriv (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    ContinuousOn (PhihetDeriv θ c w) (Set.Ioi (bHet c w)) := by
  unfold PhihetDeriv
  exact continuousOn_finsetSum _ fun i _ => continuousOn_const.mul (continuousOn_ghetDeriv hc hw i)

theorem continuousOn_PsihetDeriv (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    ContinuousOn (PsihetDeriv c w) (Set.Ioi (bHet c w)) := by
  unfold PsihetDeriv
  exact continuousOn_finsetSum _ fun i _ => continuousOn_const.mul (continuousOn_ghetDeriv hc hw i)

/-- `γ(z) = -1/s(z)` is strictly increasing on `(bHet, ∞)`. -/
theorem gamHet_lt_gamHet (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z₁ z₂ : ℝ}
    (h₁ : bHet c w < z₁) (h₁₂ : z₁ < z₂) : gamHet c w z₁ < gamHet c w z₂ := by
  have h₂ : bHet c w < z₂ := lt_trans h₁ h₁₂
  have hs₂ := sPhys_neg hc hw h₂
  have hlt := sPhys_strictMonoOn hc hw h₁ h₂ h₁₂
  unfold gamHet
  rw [neg_div, neg_div, neg_lt_neg_iff]
  exact one_div_lt_one_div_of_neg_of_lt hs₂ hlt

/-- `1 + F < 0` on `(bHet, ρ)`: `γ(z) < γ₁` and the paper's secular function is strictly
increasing above `max_i w_i²`. -/
theorem one_add_F_neg (hc : ∀ i, 0 < c i) (h4 : Assumption4 θ c w) {z : ℝ}
    (hz : bHet c w < z) (hzρ : z < rhoHet θ c w) : 1 + Fhet θ c w z < 0 := by
  have hw := exists_w_ne_zero_of_root h4.1
  have hroot := gammaTop_of_exists h4.1
  rw [one_add_F_eq hc hw hz]
  have hz0 : 0 < z := lt_trans (bHet_pos hc hw) hz
  have hγW : wSqMax w < gamHet c w z := wSqMax_lt_gamHet hc hw hz
  have hγ : 0 < gamHet c w z := lt_of_le_of_lt (wSqMax_nonneg w) hγW
  have hγlt : gamHet c w z < gammaTop θ w := by
    rw [← gamHet_rhoHet hc h4]
    exact gamHet_lt_gamHet hc hw hz hzρ
  have hsec : secular θ w (gamHet c w z) < 0 := by
    have hmono := secular_strictMonoOn (exists_signal_of_root hroot) (Set.mem_Ioi.2 hγW)
      (Set.mem_Ioi.2 hroot.1) hγlt
    rwa [hroot.2] at hmono
  exact mul_neg_of_pos_of_neg (div_pos hγ hz0) hsec

/-- `1 + F > 0` on `(ρ, ∞)`. -/
theorem one_add_F_pos (hc : ∀ i, 0 < c i) (h4 : Assumption4 θ c w) {z : ℝ}
    (hzρ : rhoHet θ c w < z) : 0 < 1 + Fhet θ c w z := by
  have hw := exists_w_ne_zero_of_root h4.1
  have hroot := gammaTop_of_exists h4.1
  have hbρ := bHet_lt_rhoHet hc h4
  have hz : bHet c w < z := lt_trans hbρ hzρ
  rw [one_add_F_eq hc hw hz]
  have hz0 : 0 < z := lt_trans (bHet_pos hc hw) hz
  have hγW : wSqMax w < gamHet c w z := wSqMax_lt_gamHet hc hw hz
  have hγ : 0 < gamHet c w z := lt_of_le_of_lt (wSqMax_nonneg w) hγW
  have hγgt : gammaTop θ w < gamHet c w z := by
    rw [← gamHet_rhoHet hc h4]
    exact gamHet_lt_gamHet hc hw hbρ hzρ
  have hsec : 0 < secular θ w (gamHet c w z) := by
    have hmono := secular_strictMonoOn (exists_signal_of_root hroot) (Set.mem_Ioi.2 hroot.1)
      (Set.mem_Ioi.2 hγW) hγgt
    rwa [hroot.2] at hmono
  exact mul_pos (div_pos hγ hz0) hsec

theorem LwDen_pos (h : ∃ g, IsGammaTop θ w g) : 0 < LwDen θ w := by
  obtain ⟨j, hj⟩ := exists_signal_of_root (gammaTop_of_exists h)
  unfold LwDen
  refine mul_pos (gammaTop_pos h) (Finset.sum_pos' (fun i _ => ?_) ⟨j, Finset.mem_univ _, ?_⟩)
  · exact div_nonneg (by positivity) (pow_nonneg (gammaTop_sub_pos h i).le 2)
  · have h1 : 0 < θ j ^ 2 * w j ^ 2 := by
      rw [show θ j ^ 2 * w j ^ 2 = (θ j * w j) ^ 2 by ring]
      exact lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 hj))
    exact div_pos h1 (pow_pos (gammaTop_sub_pos h j) 2)

/-- `F'(ρ) > 0`: `ρ F'(ρ) = LwDen/η₁` with every factor positive. -/
theorem FhetDeriv_rhoHet_pos (hc : ∀ i, 0 < c i) (h4 : Assumption4 θ c w) :
    0 < FhetDeriv θ c w (rhoHet θ c w) := by
  have hw := exists_w_ne_zero_of_root h4.1
  have hρ : 0 < rhoHet θ c w := lt_trans (bHet_pos hc hw) (bHet_lt_rhoHet hc h4)
  have h1 : 0 < rhoHet θ c w * FhetDeriv θ c w (rhoHet θ c w) := by
    rw [rhoHet_mul_FhetDeriv hc h4]
    exact div_pos (LwDen_pos h4.1) ((eta1_pos_iff θ c w).mpr h4.2)
  exact (mul_pos_iff_of_pos_left hρ).mp h1

/-- **The concrete instance of `HetScalarFacts`** (choice 1 of the header). For every edge
bound `b` with `bHet ≤ b < rhoHet`, the closed forms of `MPhet.lean` satisfy the seven
scalar statements of the analytic core. -/
theorem hetScalarFacts_of_assumption4 (hc : ∀ i, 0 < c i) (h4 : Assumption4 θ c w)
    {b : ℝ} (hb : bHet c w ≤ b) (hbρ : b < rhoHet θ c w) :
    MultiTableModel.HetScalarFacts θ c w b (rhoHet θ c w) (Phihet θ c w) (Psihet c w)
      (PhihetDeriv θ c w) (PsihetDeriv c w) where
  lt := hbρ
  neg := fun z hz hzρ => by
    have := one_add_F_neg hc h4 (lt_of_le_of_lt hb hz) hzρ
    unfold Fhet at this
    linarith
  pos := fun z hz => by
    have := one_add_F_pos hc h4 hz
    unfold Fhet at this
    linarith
  cont := by
    have hw := exists_w_ne_zero_of_root h4.1
    exact ((continuousOn_Phihet hc hw).add (continuousOn_Psihet hc hw)).mono
      (Set.Ioi_subset_Ioi hb)
  cont2 := by
    have hw := exists_w_ne_zero_of_root h4.1
    exact ((continuousOn_PhihetDeriv hc hw).add (continuousOn_PsihetDeriv hc hw)).mono
      (Set.Ioi_subset_Ioi hb)
  den_pos := FhetDeriv_rhoHet_pos hc h4
  overlap := overlap_identity_het hc h4

end MPhet

/-! ### The choice of `δ` for the overlap bracket -/

namespace R5het

/-- **How small `δ` must be** (R5's `exists_delta`, denominator only). The bracket of the
overlap is `1/((ρ+δ) D(ρ-δ)) ≤ overlap ≤ 1/((ρ-δ) D(ρ+δ))` with `D = Φ' + Ψ'`; both ends
tend to `1/(ρ D(ρ)) = L` as `δ ↓ 0`, and the proof fixes `δ` before it takes the limit in
`N`. -/
theorem exists_delta {b rho L : ℝ} {Phi2 Psi2 : ℝ → ℝ} (hbρ : b < rho) (hρ : 0 < rho)
    (hcont : ContinuousOn (fun z => Phi2 z + Psi2 z) (Set.Ioi b))
    (hD : 0 < Phi2 rho + Psi2 rho) (hL : 1 / (rho * (Phi2 rho + Psi2 rho)) = L)
    {ε : ℝ} (hε : 0 < ε) :
    ∃ δ : ℝ, 0 < δ ∧ δ < rho - b ∧ δ < rho ∧
      0 < Phi2 (rho - δ) + Psi2 (rho - δ) ∧ 0 < Phi2 (rho + δ) + Psi2 (rho + δ) ∧
      |1 / ((rho + δ) * (Phi2 (rho - δ) + Psi2 (rho - δ))) - L| < ε ∧
      |1 / ((rho - δ) * (Phi2 (rho + δ) + Psi2 (rho + δ))) - L| < ε := by
  have hwin : Set.Ioo (0 : ℝ) (min (rho - b) rho) ∈ 𝓝[>] (0 : ℝ) :=
    Ioo_mem_nhdsGT (lt_min (by linarith) hρ)
  have hev0 : ∀ᶠ δ in 𝓝[>] (0 : ℝ), 0 < δ ∧ δ < min (rho - b) rho := hwin
  have hev : ∀ᶠ δ in 𝓝[>] (0 : ℝ), 0 < δ ∧ δ < rho - b ∧ δ < rho :=
    hev0.mono fun δ h => ⟨h.1, lt_of_lt_of_le h.2 (min_le_left _ _),
      lt_of_lt_of_le h.2 (min_le_right _ _)⟩
  have hsub : Tendsto (fun δ : ℝ => rho - δ) (𝓝[>] (0 : ℝ)) (𝓝 rho) := by
    have h : Tendsto (fun δ : ℝ => rho - δ) (𝓝 (0 : ℝ)) (𝓝 (rho - 0)) :=
      tendsto_const_nhds.sub tendsto_id
    simpa using h.mono_left nhdsWithin_le_nhds
  have hadd : Tendsto (fun δ : ℝ => rho + δ) (𝓝[>] (0 : ℝ)) (𝓝 rho) := by
    have h : Tendsto (fun δ : ℝ => rho + δ) (𝓝 (0 : ℝ)) (𝓝 (rho + 0)) :=
      tendsto_const_nhds.add tendsto_id
    simpa using h.mono_left nhdsWithin_le_nhds
  have tlo : Tendsto (fun δ => Phi2 (rho - δ) + Psi2 (rho - δ)) (𝓝[>] (0 : ℝ))
      (𝓝 (Phi2 rho + Psi2 rho)) := by
    refine Tendsto.comp (hcont rho (Set.mem_Ioi.2 hbρ)) ?_
    rw [tendsto_nhdsWithin_iff]
    exact ⟨hsub, hev.mono fun δ h => Set.mem_Ioi.2 (by linarith [h.2.1])⟩
  have thi : Tendsto (fun δ => Phi2 (rho + δ) + Psi2 (rho + δ)) (𝓝[>] (0 : ℝ))
      (𝓝 (Phi2 rho + Psi2 rho)) := by
    refine Tendsto.comp (hcont rho (Set.mem_Ioi.2 hbρ)) ?_
    rw [tendsto_nhdsWithin_iff]
    exact ⟨hadd, hev.mono fun δ h => Set.mem_Ioi.2 (by linarith [h.1])⟩
  have hden0 : rho * (Phi2 rho + Psi2 rho) ≠ 0 := (mul_pos hρ hD).ne'
  have t1 : Tendsto (fun δ => 1 / ((rho + δ) * (Phi2 (rho - δ) + Psi2 (rho - δ))))
      (𝓝[>] (0 : ℝ)) (𝓝 L) := by
    have h := (tendsto_const_nhds (x := (1 : ℝ))).div (hadd.mul tlo) hden0
    rwa [hL] at h
  have t2 : Tendsto (fun δ => 1 / ((rho - δ) * (Phi2 (rho + δ) + Psi2 (rho + δ))))
      (𝓝[>] (0 : ℝ)) (𝓝 L) := by
    have h := (tendsto_const_nhds (x := (1 : ℝ))).div (hsub.mul thi) hden0
    rwa [hL] at h
  have e1 : ∀ᶠ δ in 𝓝[>] (0 : ℝ), 0 < Phi2 (rho - δ) + Psi2 (rho - δ) :=
    tlo.eventually_const_lt hD
  have e2 : ∀ᶠ δ in 𝓝[>] (0 : ℝ), 0 < Phi2 (rho + δ) + Psi2 (rho + δ) :=
    thi.eventually_const_lt hD
  have e3 := t1.eventually (Metric.ball_mem_nhds L hε)
  have e4 := t2.eventually (Metric.ball_mem_nhds L hε)
  obtain ⟨δ, ⟨hδ0, hδb, hδρ⟩, h1, h2, h3, h4⟩ := (hev.and (e1.and (e2.and (e3.and e4)))).exists
  refine ⟨δ, hδ0, hδb, hδρ, h1, h2, ?_, ?_⟩
  · simpa only [Metric.mem_ball, Real.dist_eq] using h3
  · simpa only [Metric.mem_ball, Real.dist_eq] using h4

end R5het

/-! ### The model layer -/

namespace MultiTableModel

open R4

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]

/-! #### Measurability of the edge event (audit item F2) -/

/-- `W₀' = Y Yᵀ` with `Y = Σ^{1/2} E⊥`. -/
theorem W0het_eq_mul_transpose (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) :
    m.W0het w N ω
      = (m.SigmaHalf w N * m.EperpHet N ω) * (m.SigmaHalf w N * m.EperpHet N ω)ᵀ := by
  rw [W0het, Matrix.transpose_mul, transpose_SigmaHalf, Matrix.mul_assoc, Matrix.mul_assoc]

theorem measurable_EperpHet (m : MultiTableModel μ M n d) (N : ℕ) :
    Measurable (m.EperpHet N) := by
  have hZ : ∀ (r : Fin (∑ i, n i N)) (j : Fin (d N)), Measurable fun ω => m.stack.Z N ω r j :=
    fun r j => (measurable_pi_apply j).comp ((measurable_pi_apply r).comp (m.stack.hZ N))
  have hE : ∀ (r : Fin (∑ i, n i N)) (j : Fin (d N)), Measurable fun ω => m.stack.E N ω r j := by
    intro r j
    have h : (fun ω => m.stack.E N ω r j)
        = fun ω => (Real.sqrt (d N))⁻¹ * m.stack.Z N ω r j := rfl
    rw [h]
    exact (hZ r j).const_mul _
  have he : ∀ r : Fin (∑ i, n i N), Measurable fun ω => m.eHet N ω r := by
    intro r
    have h : (fun ω => m.eHet N ω r)
        = fun ω => ∑ j, m.stack.E N ω r j * WithLp.ofLp (m.stack.v N) j := rfl
    rw [h]
    exact Finset.measurable_sum _ fun j _ => (hE r j).mul_const _
  refine measurable_pi_lambda _ fun r => measurable_pi_lambda _ fun j => ?_
  have h : (fun ω => m.EperpHet N ω r j)
      = fun ω => m.stack.E N ω r j - m.eHet N ω r * WithLp.ofLp (m.stack.v N) j := rfl
  rw [h]
  exact (hE r j).sub ((he r).mul_const _)

theorem measurable_Yhet_transpose (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) :
    Measurable fun ω => (m.SigmaHalf w N * m.EperpHet N ω)ᵀ := by
  refine measurable_pi_lambda _ fun j => measurable_pi_lambda _ fun r => ?_
  have h : (fun ω => (m.SigmaHalf w N * m.EperpHet N ω)ᵀ j r)
      = fun ω => ∑ k, m.SigmaHalf w N r k * m.EperpHet N ω k j := rfl
  rw [h]
  exact Finset.measurable_sum _ fun k _ =>
    ((measurable_pi_apply j).comp ((measurable_pi_apply k).comp
      (m.measurable_EperpHet N))).const_mul _

/-- `λ_max(W₀') = gramLamMax Yᵀ`: `W₀' = Y Yᵀ = (Yᵀ)ᵀ Yᵀ`. -/
theorem lamMax_W0het_eq_gramLamMax (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) :
    lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω)
      = gramLamMax (m.SigmaHalf w N * m.EperpHet N ω)ᵀ := by
  have h : m.W0het w N ω = ((m.SigmaHalf w N * m.EperpHet N ω)ᵀ)ᵀ
      * (m.SigmaHalf w N * m.EperpHet N ω)ᵀ := by
    rw [W0het_eq_mul_transpose, Matrix.transpose_transpose]
  exact lamMax_congr h _ _

/-- The edge event is measurable (audit item F2). -/
theorem measurableSet_lamMax_W0het_le (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ)
    (r : ℝ) :
    MeasurableSet {ω | lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω) ≤ r} := by
  have h : {ω | lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω) ≤ r}
      = (fun ω => gramLamMax (m.SigmaHalf w N * m.EperpHet N ω)ᵀ) ⁻¹' Set.Iic r := by
    ext ω
    simp only [Set.mem_ofPred_eq, Set.mem_preimage, Set.mem_Iic, lamMax_W0het_eq_gramLamMax]
  rw [h]
  exact (measurable_gramLamMax.comp (m.measurable_Yhet_transpose w N)) measurableSet_Iic

theorem lamMax_W0het_nonneg (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    0 ≤ lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω) :=
  Het.lamMax_nonneg_of_posSemidef _ (m.posSemidef_W0het w N ω) (m.stack.hn N)

/-! #### Localization at the secular root (R5 step 4, on the `n` side) -/

/-- **One `ω`.** If `λ_max(W₀') < zlo`, the secular function of `W₀' + q qᵀ = X_w X_wᵀ` is
negative at `zlo` and positive at `zhi`, then `q ≠ 0`, the top eigenvalue of `X_wᵀ X_w`
(equal to that of `X_w X_wᵀ`) is the secular root inside `(zlo, zhi)`, and it is simple on
the `d` side too. -/
theorem localize_het (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N)
    {zlo zhi : ℝ} (h₀ : lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω) < zlo)
    (hz : zlo < zhi) (hslo : secular (m.W0het w N ω) (m.qHet w N ω) zlo < 0)
    (hshi : 0 < secular (m.W0het w N ω) (m.qHet w N ω) zhi) :
    m.qHet w N ω ≠ 0 ∧ zlo < gramLamMax ((m.stackW w).X N ω) ∧
      gramLamMax ((m.stackW w).X N ω) < zhi ∧
      secular (m.W0het w N ω) (m.qHet w N ω) (gramLamMax ((m.stackW w).X N ω)) = 0 ∧
      TopSimple (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) := by
  have hW := m.isHermitian_W0het w N ω
  have hq : m.qHet w N ω ≠ 0 := by
    intro h0
    have h1 : secular (m.W0het w N ω) (m.qHet w N ω) zlo = 1 := by
      change 1 + qform (m.W0het w N ω) zlo (m.qHet w N ω) = 1
      rw [h0]
      simp [qform]
    linarith
  have hcont : ContinuousOn (fun t => secular (m.W0het w N ω) (m.qHet w N ω) t)
      (Set.Icc zlo zhi) := by
    intro s hs
    have hs' : lamMax (m.W0het w N ω) hW < s := lt_of_lt_of_le h₀ hs.1
    exact (continuousAt_const.add
      (hasDerivAt_qform (y := m.qHet w N ω) hW hs').continuousAt).continuousWithinAt
  obtain ⟨lam, hlmem, hlz⟩ :=
    intermediate_value_Ioo hz.le hcont (Set.mem_Ioo.2 ⟨hslo, hshi⟩)
  have hlamlt : lamMax (m.W0het w N ω) hW < lam := lt_trans h₀ hlmem.1
  have hA' : (m.W0het w N ω + Matrix.vecMulVec (m.qHet w N ω) (m.qHet w N ω)).IsHermitian := by
    rw [← m.gram_eq_het w N ω]; exact Het.isHermitian_mul_transpose_self _
  have hgl : gramLamMax ((m.stackW w).X N ω)
      = lamMax (m.W0het w N ω + Matrix.vecMulVec (m.qHet w N ω) (m.qHet w N ω)) hA' := by
    rw [Het.lamMax_gram_comm (m.stack.hn N) (m.stack.hd N)]
    exact lamMax_congr (m.gram_eq_het w N ω) _ hA'
  have hlamMax : gramLamMax ((m.stackW w).X N ω) = lam := by
    rw [hgl, lamMax_eq hW hq hlamlt hlz hA']
  have hsimpleN : TopSimple ((m.stackW w).X N ω * ((m.stackW w).X N ω)ᵀ)
      (Het.isHermitian_mul_transpose_self _) :=
    topSimple_congr (m.gram_eq_het w N ω) _ hA' (topSimple hW hq hlamlt hlz hA')
  have hpos : 0 < gramLamMax ((m.stackW w).X N ω) := by
    rw [hlamMax]; exact lt_of_le_of_lt (m.lamMax_W0het_nonneg w N ω) hlamlt
  refine ⟨hq, ?_, ?_, ?_, ?_⟩
  · rw [hlamMax]; exact hlmem.1
  · rw [hlamMax]; exact hlmem.2
  · rw [hlamMax]; exact hlz
  · exact Het.topSimple_transpose_mul_self (m.stack.hn N) (m.stack.hd N) _ hpos hsimpleN

/-- **The overlap bracket, one `ω`** (R5 step 7 without the numerator). On the localization
event, `overlap = 1/(λ · Q(λ))` with `Q(z) = qform2 W₀' z q` decreasing in `z`, so
`1/(zhi · Q(zlo)) ≤ overlap ≤ 1/(zlo · Q(zhi))`. -/
theorem overlap_bracket_het (m : MultiTableModel μ M n d) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N)
    {zlo zhi : ℝ} (hzlo : 0 < zlo)
    (h₀ : lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω) < zlo) (hz : zlo < zhi)
    (hslo : secular (m.W0het w N ω) (m.qHet w N ω) zlo < 0)
    (hshi : 0 < secular (m.W0het w N ω) (m.qHet w N ω) zhi) :
    1 / (zhi * qform2 (m.W0het w N ω) zlo (m.qHet w N ω)) ≤ m.stackPerfW w N ω ∧
      m.stackPerfW w N ω ≤ 1 / (zlo * qform2 (m.W0het w N ω) zhi (m.qHet w N ω)) := by
  obtain ⟨hq, hlo, hhi, -, -⟩ := m.localize_het w N ω h₀ hz hslo hshi
  have hW := m.isHermitian_W0het w N ω
  set lam := gramLamMax ((m.stackW w).X N ω) with hlamdef
  have hlamlt : lamMax (m.W0het w N ω) hW < lam := lt_trans h₀ hlo
  rw [m.stackPerfW_eq_inv_lam_qform2 w N ω hlamlt]
  have hQpos : ∀ z, lamMax (m.W0het w N ω) hW < z →
      0 < qform2 (m.W0het w N ω) z (m.qHet w N ω) := by
    intro z hz'
    unfold qform2
    rw [dotProduct_resolv_sq_eq_norm_sq hW]
    exact pow_pos (norm_pos_iff.mpr (resolv_toLp_ne_zero hW hz' hq)) 2
  have hQlo := hQpos zlo h₀
  have hQhi := hQpos zhi (lt_trans hlamlt hhi)
  have hQlam := hQpos lam hlamlt
  have h1 : qform2 (m.W0het w N ω) lam (m.qHet w N ω)
      ≤ qform2 (m.W0het w N ω) zlo (m.qHet w N ω) :=
    qform2_antitoneOn hW _ (Set.mem_Ioi.2 h₀) (Set.mem_Ioi.2 hlamlt) hlo.le
  have h2 : qform2 (m.W0het w N ω) zhi (m.qHet w N ω)
      ≤ qform2 (m.W0het w N ω) lam (m.qHet w N ω) :=
    qform2_antitoneOn hW _ (Set.mem_Ioi.2 hlamlt) (Set.mem_Ioi.2 (lt_trans hlamlt hhi)) hhi.le
  have hlampos : 0 < lam := lt_trans hzlo hlo
  have hprod1 : zlo * qform2 (m.W0het w N ω) zhi (m.qHet w N ω)
      ≤ lam * qform2 (m.W0het w N ω) lam (m.qHet w N ω) :=
    mul_le_mul hlo.le h2 hQhi.le hlampos.le
  have hprod2 : lam * qform2 (m.W0het w N ω) lam (m.qHet w N ω)
      ≤ zhi * qform2 (m.W0het w N ω) zlo (m.qHet w N ω) :=
    mul_le_mul hhi.le h1 hQlam.le (by linarith)
  exact ⟨one_div_le_one_div_of_le (mul_pos hlampos hQlam) hprod2,
    one_div_le_one_div_of_le (mul_pos hzlo hQhi) hprod1⟩

/-! #### Limits of the random forms (R5 steps 2 and 3) -/

variable [∀ N, IsProbabilityMeasure (μ N)] {m : MultiTableModel μ M n d} {w c : Fin M → ℝ}
  {b rho : ℝ} {Phi Psi Phi2 Psi2 : ℝ → ℝ}

set_option linter.unusedSectionVars false

/-- The random secular function converges at every fixed `z > b`, to `1 + Φ(z) + Ψ(z)`. -/
theorem tendstoInProb_secular_het (H : m.ResolventLimitsHet w b Phi Psi Phi2 Psi2) {z : ℝ}
    (hz : b < z) :
    TendstoInProb μ (fun N ω => secular (m.W0het w N ω) (m.qHet w N ω) z)
      (1 + Phi z + Psi z) := by
  have hfun : ∀ (N : ℕ) (ω : Ω N), secular (m.W0het w N ω) (m.qHet w N ω) z
      = 1 + qform (m.W0het w N ω) z (m.u0Het w N)
          + 2 * cform (m.W0het w N ω) z (m.u0Het w N) (m.SigmaHalf w N *ᵥ m.eHet N ω)
          + qform (m.W0het w N ω) z (m.SigmaHalf w N *ᵥ m.eHet N ω) := by
    intro N ω
    change 1 + qform (m.W0het w N ω) z (m.qHet w N ω) = _
    rw [show m.qHet w N ω = m.u0Het w N + m.SigmaHalf w N *ᵥ m.eHet N ω from rfl,
      qform_add (m.isHermitian_W0het w N ω)]
    ring
  have h := (((TendstoInProb.const μ (1 : ℝ)).add (H.uu z hz)).add
    ((H.ue z hz).const_mul 2)).add (H.ee z hz)
  have hlim : (1 : ℝ) + Phi z + 2 * 0 + Psi z = 1 + Phi z + Psi z := by ring
  rw [hlim] at h
  simpa only [hfun] using h

/-- The denominator form at `q` converges at every fixed `z > b`, to `Φ'(z) + Ψ'(z)`. -/
theorem tendstoInProb_qform2_qHet (H : m.ResolventLimitsHet w b Phi Psi Phi2 Psi2) {z : ℝ}
    (hz : b < z) :
    TendstoInProb μ (fun N ω => qform2 (m.W0het w N ω) z (m.qHet w N ω))
      (Phi2 z + Psi2 z) := by
  have hfun : ∀ (N : ℕ) (ω : Ω N), qform2 (m.W0het w N ω) z (m.qHet w N ω)
      = qform2 (m.W0het w N ω) z (m.u0Het w N)
          + 2 * cform2 (m.W0het w N ω) z (m.u0Het w N) (m.SigmaHalf w N *ᵥ m.eHet N ω)
          + qform2 (m.W0het w N ω) z (m.SigmaHalf w N *ᵥ m.eHet N ω) := by
    intro N ω
    rw [show m.qHet w N ω = m.u0Het w N + m.SigmaHalf w N *ᵥ m.eHet N ω from rfl,
      qform2_add (m.isHermitian_W0het w N ω)]
  have h := ((H.uu2 z hz).add ((H.ue2 z hz).const_mul 2)).add (H.ee2 z hz)
  have hlim : Phi2 z + 2 * (0 : ℝ) + Psi2 z = Phi2 z + Psi2 z := by ring
  rw [hlim] at h
  simpa only [hfun] using h

/-- **The localization bad set.** Its measure tends to `0`, and off it the hypotheses of
`localize_het` hold at `zlo = ρ - δ` and `zhi = ρ + δ`. The edge field enters here only. -/
theorem exists_localization_bad_het (H : m.ResolventLimitsHet w b Phi Psi Phi2 Psi2)
    (hsc : HetScalarFacts (fun i => (m.tbl i).θ) c w b rho Phi Psi Phi2 Psi2)
    {δ : ℝ} (hδ : 0 < δ) (hδ' : δ < rho - b) :
    ∃ B : ∀ N, Set (Ω N), Tendsto (fun N => μ N (B N)) atTop (𝓝 0) ∧
      ∀ (N : ℕ) (ω : Ω N), ω ∉ B N →
        lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω) < rho - δ ∧
        secular (m.W0het w N ω) (m.qHet w N ω) (rho - δ) < 0 ∧
        0 < secular (m.W0het w N ω) (m.qHet w N ω) (rho + δ) := by
  have hbρ := hsc.lt
  have hlo : b < rho - δ := by linarith
  have hhi : b < rho + δ := by linarith
  have hFlo : 1 + Phi (rho - δ) + Psi (rho - δ) < 0 := hsc.neg _ hlo (by linarith)
  have hFhi : 0 < 1 + Phi (rho + δ) + Psi (rho + δ) := hsc.pos _ (by linarith)
  have hε₀ : (0 : ℝ) < (rho - δ - b) / 2 := by linarith
  have hB1 : Tendsto (fun N => μ N
      {ω | lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω)
        ≤ b + (rho - δ - b) / 2}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_lamMax_W0het_le w N _).nullMeasurableSet) (H.edge _ hε₀)
  have hB2 := tendstoInProb_secular_het H hlo
    (-(1 + Phi (rho - δ) + Psi (rho - δ)) / 2) (by linarith)
  have hB3 := tendstoInProb_secular_het H hhi
    ((1 + Phi (rho + δ) + Psi (rho + δ)) / 2) (by linarith)
  refine ⟨_, tendsto_measure_zero_union hB1 (tendsto_measure_zero_union hB2 hB3), ?_⟩
  intro N ω hω
  simp only [Set.mem_union, not_or] at hω
  obtain ⟨g1, g2, g3⟩ := hω
  have hlam : lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω)
      ≤ b + (rho - δ - b) / 2 := by
    by_contra hx; exact g1 (Set.mem_compl hx)
  refine ⟨by linarith, ?_, ?_⟩
  · by_contra hx
    rw [not_lt] at hx
    refine g2 ?_
    change -(1 + Phi (rho - δ) + Psi (rho - δ)) / 2
      ≤ |secular (m.W0het w N ω) (m.qHet w N ω) (rho - δ)
        - (1 + Phi (rho - δ) + Psi (rho - δ))|
    rw [abs_of_nonneg (by linarith)]
    linarith
  · by_contra hx
    rw [not_lt] at hx
    refine g3 ?_
    change (1 + Phi (rho + δ) + Psi (rho + δ)) / 2
      ≤ |secular (m.W0het w N ω) (m.qHet w N ω) (rho + δ)
        - (1 + Phi (rho + δ) + Psi (rho + δ))|
    rw [abs_of_nonpos (by linarith)]
    linarith

/-- The edge bound is nonnegative: `W₀'` is positive semidefinite, so an edge event at a
negative level is empty and could not fill the space (choice 2 of the header). -/
theorem b_nonneg_of_edge (H : m.ResolventLimitsHet w b Phi Psi Phi2 Psi2) : 0 ≤ b := by
  by_contra hb
  rw [not_le] at hb
  have h := H.edge (-b / 2) (by linarith)
  have hempty : ∀ N, {ω | lamMax (m.W0het w N ω) (m.isHermitian_W0het w N ω) ≤ b + -b / 2}
      = (∅ : Set (Ω N)) := by
    intro N
    ext ω
    simp only [Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false, not_le]
    linarith [m.lamMax_W0het_nonneg w N ω]
  simp only [hempty, measure_empty] at h
  exact one_ne_zero (tendsto_nhds_unique h tendsto_const_nhds)

/-! #### The three conclusions of item H4 -/

/-- **H4 conclusion 1.** `λ_max(X_wᵀ X_w) → ρ` in probability. -/
theorem lamMax_tendstoInProb_het (H : m.ResolventLimitsHet w b Phi Psi Phi2 Psi2)
    (hsc : HetScalarFacts (fun i => (m.tbl i).θ) c w b rho Phi Psi Phi2 Psi2) :
    TendstoInProb μ (fun N ω => gramLamMax ((m.stackW w).X N ω)) rho := by
  have hbρ := hsc.lt
  intro ε hε
  obtain ⟨δ, hδ, hδ', hδε⟩ : ∃ δ : ℝ, 0 < δ ∧ δ < rho - b ∧ δ ≤ ε :=
    ⟨min ε ((rho - b) / 2), lt_min hε (by linarith),
      lt_of_le_of_lt (min_le_right _ _) (by linarith), min_le_left _ _⟩
  obtain ⟨B, hB, hBgood⟩ := exists_localization_bad_het H hsc hδ hδ'
  refine tendsto_measure_zero_of_subset ?_ hB
  intro N ω hω
  by_contra hnb
  obtain ⟨hlam, hslo, hshi⟩ := hBgood N ω hnb
  obtain ⟨-, hlo, hhi, -, -⟩ := m.localize_het w N ω hlam (by linarith) hslo hshi
  have h2 : ε ≤ |gramLamMax ((m.stackW w).X N ω) - rho| := hω
  have h3 : |gramLamMax ((m.stackW w).X N ω) - rho| < δ := by
    rw [abs_lt]; constructor <;> linarith
  linarith

/-- **H4 conclusion 2.** The top eigenvalue of `X_wᵀ X_w` is simple with probability tending
to `1`. The almost sure version is item H10, not this file. -/
theorem tendsto_measure_topSimple_het (H : m.ResolventLimitsHet w b Phi Psi Phi2 Psi2)
    (hsc : HetScalarFacts (fun i => (m.tbl i).θ) c w b rho Phi Psi Phi2 Psi2) :
    Tendsto (fun N => μ N {ω | TopSimple (m.stackGramW w N ω)
      (m.isHermitian_stackGramW w N ω)}) atTop (𝓝 1) := by
  have hbρ := hsc.lt
  obtain ⟨δ, hδ, hδ'⟩ : ∃ δ : ℝ, 0 < δ ∧ δ < rho - b :=
    ⟨(rho - b) / 2, by linarith, by linarith⟩
  obtain ⟨B, hB, hBgood⟩ := exists_localization_bad_het H hsc hδ hδ'
  refine tendsto_measure_one_of_bad ?_ hB
  intro N ω hω
  by_contra hnb
  obtain ⟨hlam, hslo, hshi⟩ := hBgood N ω hnb
  obtain ⟨-, -, -, -, hsimple⟩ := m.localize_het w N ω hlam (by linarith) hslo hshi
  exact hω hsimple

/-- **H4 conclusion 3.** `overlap(X_w, v) → L(w)` in probability, the `align` field of
`HeteroLaw`. The `ε`-`δ` order is R5's: `δ` is fixed first (`R5het.exists_delta`), then
the limit in `N` is taken at the two points `ρ ∓ δ`. -/
theorem align_tendstoInProb_het (H : m.ResolventLimitsHet w b Phi Psi Phi2 Psi2)
    (hsc : HetScalarFacts (fun i => (m.tbl i).θ) c w b rho Phi Psi Phi2 Psi2) :
    TendstoInProb μ (fun N ω => m.stackPerfW w N ω)
      (Scalars.Lw (fun i => (m.tbl i).θ) c w) := by
  have hb0 : 0 ≤ b := b_nonneg_of_edge H
  have hρ : 0 < rho := lt_of_le_of_lt hb0 hsc.lt
  intro ε hε
  obtain ⟨δ, hδ, hδ', hδρ, hDlo, hDhi, hLo, hHi⟩ :=
    R5het.exists_delta hsc.lt hρ hsc.cont2 hsc.den_pos hsc.overlap (half_pos hε)
  have hlo : b < rho - δ := by linarith
  have hhi : b < rho + δ := by linarith
  obtain ⟨B, hB, hBgood⟩ := exists_localization_bad_het H hsc hδ hδ'
  have hne1 : (rho + δ) * (Phi2 (rho - δ) + Psi2 (rho - δ)) ≠ 0 :=
    (mul_pos (by linarith) hDlo).ne'
  have hne2 : (rho - δ) * (Phi2 (rho + δ) + Psi2 (rho + δ)) ≠ 0 :=
    (mul_pos (by linarith) hDhi).ne'
  have hB2 := ((TendstoInProb.const μ (1 : ℝ)).div
    ((TendstoInProb.const μ (rho + δ)).mul (tendstoInProb_qform2_qHet H hlo)) hne1)
    (ε / 2) (half_pos hε)
  have hB3 := ((TendstoInProb.const μ (1 : ℝ)).div
    ((TendstoInProb.const μ (rho - δ)).mul (tendstoInProb_qform2_qHet H hhi)) hne2)
    (ε / 2) (half_pos hε)
  refine tendsto_measure_zero_of_subset ?_
    (tendsto_measure_zero_union hB (tendsto_measure_zero_union hB2 hB3))
  intro N ω hω
  by_contra hcon
  simp only [Set.mem_union, not_or] at hcon
  obtain ⟨g0, g2, g3⟩ := hcon
  obtain ⟨hlam, hslo, hshi⟩ := hBgood N ω g0
  have g2' : |1 / ((rho + δ) * qform2 (m.W0het w N ω) (rho - δ) (m.qHet w N ω))
      - 1 / ((rho + δ) * (Phi2 (rho - δ) + Psi2 (rho - δ)))| < ε / 2 := not_le.mp g2
  have g3' : |1 / ((rho - δ) * qform2 (m.W0het w N ω) (rho + δ) (m.qHet w N ω))
      - 1 / ((rho - δ) * (Phi2 (rho + δ) + Psi2 (rho + δ)))| < ε / 2 := not_le.mp g3
  obtain ⟨hlow, hhigh⟩ :=
    m.overlap_bracket_het w N ω (by linarith) hlam (by linarith) hslo hshi
  have h2 := abs_lt.mp g2'
  have h3 := abs_lt.mp g3'
  have h8 := abs_lt.mp hLo
  have h9 := abs_lt.mp hHi
  have hfin : |m.stackPerfW w N ω - Scalars.Lw (fun i => (m.tbl i).θ) c w| < ε := by
    rw [abs_lt]
    constructor <;> linarith
  exact absurd (show ε ≤ |m.stackPerfW w N ω - Scalars.Lw (fun i => (m.tbl i).θ) c w|
    from hω) (not_le.mpr hfin)

/-! #### The concrete branch (the H4 row of `notes/archive/plan_heterolaw_A.md` section 4)

`b` is any edge bound with `bHet ≤ b < rhoHet`: `b = bHet` (hypothesis `HeteroEdge`) or
`b = bSF` with the margin `bSF < rhoHet` (item H8, `bHet_le_bSF`). -/

/-- **H4 conclusion 1, concrete.** -/
theorem lamMax_tendstoInProb_het_of_assumption4 (hc : ∀ i, 0 < c i)
    (h4 : Scalars.Assumption4 (fun i => (m.tbl i).θ) c w) {b : ℝ}
    (hb : MPhet.bHet c w ≤ b) (hbρ : b < MPhet.rhoHet (fun i => (m.tbl i).θ) c w)
    (H : m.ResolventLimitsHet w b (MPhet.Phihet (fun i => (m.tbl i).θ) c w)
      (MPhet.Psihet c w) (MPhet.PhihetDeriv (fun i => (m.tbl i).θ) c w)
      (MPhet.PsihetDeriv c w)) :
    TendstoInProb μ (fun N ω => gramLamMax ((m.stackW w).X N ω))
      (MPhet.rhoHet (fun i => (m.tbl i).θ) c w) :=
  lamMax_tendstoInProb_het H (MPhet.hetScalarFacts_of_assumption4 hc h4 hb hbρ)

/-- **H4 conclusion 2, concrete.** -/
theorem tendsto_measure_topSimple_het_of_assumption4 (hc : ∀ i, 0 < c i)
    (h4 : Scalars.Assumption4 (fun i => (m.tbl i).θ) c w) {b : ℝ}
    (hb : MPhet.bHet c w ≤ b) (hbρ : b < MPhet.rhoHet (fun i => (m.tbl i).θ) c w)
    (H : m.ResolventLimitsHet w b (MPhet.Phihet (fun i => (m.tbl i).θ) c w)
      (MPhet.Psihet c w) (MPhet.PhihetDeriv (fun i => (m.tbl i).θ) c w)
      (MPhet.PsihetDeriv c w)) :
    Tendsto (fun N => μ N {ω | TopSimple (m.stackGramW w N ω)
      (m.isHermitian_stackGramW w N ω)}) atTop (𝓝 1) :=
  tendsto_measure_topSimple_het H (MPhet.hetScalarFacts_of_assumption4 hc h4 hb hbρ)

/-- **H4 conclusion 3, concrete.** -/
theorem align_tendstoInProb_het_of_assumption4 (hc : ∀ i, 0 < c i)
    (h4 : Scalars.Assumption4 (fun i => (m.tbl i).θ) c w) {b : ℝ}
    (hb : MPhet.bHet c w ≤ b) (hbρ : b < MPhet.rhoHet (fun i => (m.tbl i).θ) c w)
    (H : m.ResolventLimitsHet w b (MPhet.Phihet (fun i => (m.tbl i).θ) c w)
      (MPhet.Psihet c w) (MPhet.PhihetDeriv (fun i => (m.tbl i).θ) c w)
      (MPhet.PsihetDeriv c w)) :
    TendstoInProb μ (fun N ω => m.stackPerfW w N ω)
      (Scalars.Lw (fun i => (m.tbl i).θ) c w) :=
  align_tendstoInProb_het H (MPhet.hetScalarFacts_of_assumption4 hc h4 hb hbρ)

end MultiTableModel
end StackedSVD
