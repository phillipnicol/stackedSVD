/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.SingleWeight.Scalars
import StackedSVD.RMT.Het.MPhet

/-!
# Unit G0: the matrix secular layer of `prop:gen_rank_stacksvd_singleweight`

Track G, unit G0 of `notes/archive/trackG_plan.md` (2026-09-05). Namespace
`StackedSVD.SingleWeight`. Nothing here is random: every declaration is an identity between
finitely many real matrices and real numbers.

This file mirrors, at general rank `r`, the scalar layer of `RMT/Het/MPhet.lean`:

| here | rank-one mirror in `RMT/Het/MPhet.lean` |
|---|---|
| `swFmat`, `one_add_swFmat_eq` | `MPhet.Fhet`, `MPhet.one_add_F_eq` (`:712`) |
| `swRho`, `bHet_lt_swRho` | `MPhet.rhoHet_eq_zfun` (`:604`), `MPhet.bHet_lt_rhoHet` (`:673`) |
| `phi_neg_inv_eq` | `MPhet.phi_neg_inv_gammaTop` (`:620`) |
| `swFmat2`, `qform_swFmat2_eq_swNu` | `MPhet.FhetDeriv`, `MPhet.rhoHet_mul_FhetDeriv` (`:788`) |
| `one_div_swRho_mul_swNu` | `MPhet.overlap_identity_het` (`:842`) |

## Content

1. `swFmat`, `swFmat2`: the matrix valued secular function `I + F(z)` of the weighted stack
   and its first derivative in `z`, with the signal term carried by `sigMat`.
2. `swRho c w γ = MPhet.zfun c w (-1/γ)`: the outlier attached to a secular root `γ`.
3. `one_add_swFmat_eq`: `I + F(z) = (γ(z)/z) (I - secMat θ R w γ(z))` for every `z` above
   the bulk edge. This is the matrix form of `MPhet.one_add_F_eq`.
4. `bHet_lt_swRho`, `gamHet_swRho`, `swRho_lt_swRho`: the outlier lies above the edge, the
   physical branch takes it back to `γ`, and the map `γ ↦ swRho c w γ` is strictly
   increasing.
5. `swFmat_mulVec_eigvec`: `F(ρ) z = -z` for a unit eigenvector `z` of `secMat` at `γ`.
6. `swNu`, `qform_swFmat2_eq_swNu`, `swNu_pos`, `one_div_swRho_mul_swNu`: the scalar
   `ν = zᵀ F'(ρ) z` in closed form, its positivity, and `1/(ρ ν) = swTerm`.

## Hypotheses

`hc : ∀ i, 0 < c i`, `hw : ∃ i, w i ≠ 0`, `hγ : Scalars.wSqMax w < γ`,
`hη : ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 < 1` (the detectability threshold of
`main_paper.tex:2106`), `hznorm : ‖z‖ = 1` and `hz : secMat θ R w γ *ᵥ z = z`.
`hw` is not implied by `hγ` and `hη`: at `w = 0` the junk value of `MPhet.sPhys` breaks
`gamHet_swRho`. `exists_w_ne_zero_of_isSecularRoot` supplies `hw` from `IsSecularRoot`.

Numeric check of every identity: `notes/archive/trackG_specs/check_trackG.py`, seed 20260905, checks
1, 2, 3 and 3b, errors 1.3e-16 to 2.0e-10.
-/

open Filter Topology
open scoped Matrix

namespace StackedSVD

namespace SingleWeight

variable {M r : ℕ} {rk : Fin M → ℕ}

/-! ### 0. Elementary positivity above `max_i w_i²` -/

/-- Above `max_i w_i²` every `γ - w_i²` is positive. -/
theorem gam_sub_pos {w : Fin M → ℝ} {γ : ℝ} (hγ : Scalars.wSqMax w < γ) (i : Fin M) :
    0 < γ - w i ^ 2 := by
  have := Scalars.le_wSqMax w i
  linarith

/-- Above `max_i w_i²` the root itself is positive. -/
theorem gam_pos {w : Fin M → ℝ} {γ : ℝ} (hγ : Scalars.wSqMax w < γ) : 0 < γ :=
  lt_of_le_of_lt (Scalars.wSqMax_nonneg w) hγ

/-- A secular root forces some weight to be nonzero: at `w = 0` the matrix `secMat` is `0`
and `det (1 - 0) = 1 ≠ 0`. This is how a consumer supplies `hw`. -/
theorem exists_w_ne_zero_of_isSecularRoot {θ : (i : Fin M) → Fin (rk i) → ℝ}
    {R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ} {w : Fin M → ℝ} {γ : ℝ}
    (h : IsSecularRoot θ R w γ) : ∃ i, w i ≠ 0 := by
  by_contra hcon
  have hcon' : ∀ i, w i = 0 := fun i => by
    by_contra hi
    exact hcon ⟨i, hi⟩
  have h0 : secMat θ R w γ = 0 := by
    simp only [secMat]
    refine Finset.sum_eq_zero fun i _ => ?_
    rw [hcon' i]
    simp
  have hdet := h.2
  rw [h0, sub_zero, Matrix.det_one] at hdet
  exact one_ne_zero hdet

/-! ### 1. The two matrix valued secular functions -/

/-- `I + F(z)` without the `I`: the matrix `F(z) = ∑_i w_i² g_i(z) S_i + Ψ(z) I` of the
weighted stack, with `S_i = R_i Θ_i² R_iᵀ`. At `r = 1` and `S_i = θ_i²` it is
`MPhet.Fhet θ c w z`. -/
noncomputable def swFmat (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (c w : Fin M → ℝ) (z : ℝ) :
    Matrix (Fin r) (Fin r) ℝ :=
  (∑ i, (w i ^ 2 * MPhet.ghet c w i z) • sigMat θ R i)
    + MPhet.Psihet c w z • (1 : Matrix (Fin r) (Fin r) ℝ)

/-- `F'(z)`, the derivative of `swFmat` in `z`, in the closed forms of `MPhet.ghetDeriv` and
`MPhet.PsihetDeriv`. -/
noncomputable def swFmat2 (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (c w : Fin M → ℝ) (z : ℝ) :
    Matrix (Fin r) (Fin r) ℝ :=
  (∑ i, (w i ^ 2 * MPhet.ghetDeriv c w i z) • sigMat θ R i)
    + MPhet.PsihetDeriv c w z • (1 : Matrix (Fin r) (Fin r) ℝ)

/-- The outlier attached to a secular root `γ`: `ρ = zfun(-1/γ)`. It is the rank-`r` twin of
`MPhet.rhoHet`, which `MPhet.rhoHet_eq_zfun` writes in exactly this form. -/
noncomputable def swRho (c w : Fin M → ℝ) (γ : ℝ) : ℝ := MPhet.zfun c w (-1 / γ)

/-- `ν_ℓ = z_ℓᵀ F'(ρ_ℓ) z_ℓ` in closed form: `γ (zᵀ K(γ) z)/(η(γ) ρ)` with
`K = secDerivMat` and `η(γ) = 1 - ∑_i c_i w_i⁴/(γ - w_i²)²`. -/
noncomputable def swNu (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (c w : Fin M → ℝ) (γ : ℝ)
    (z : EuclideanSpace ℝ (Fin r)) : ℝ :=
  γ * (WithLp.ofLp z ⬝ᵥ (secDerivMat θ R w γ *ᵥ WithLp.ofLp z))
    / ((1 - ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2) * swRho c w γ)

/-! ### 2. Quadratic forms of the three `sigMat` sums -/

/-- The quadratic form of a scalar combination of the `S_i`. -/
theorem qform_sum (a : Fin M → ℝ) (S : Fin M → Matrix (Fin r) (Fin r) ℝ) (v : Fin r → ℝ) :
    v ⬝ᵥ ((∑ i, a i • S i) *ᵥ v) = ∑ i, a i * (v ⬝ᵥ (S i *ᵥ v)) := by
  rw [Matrix.sum_mulVec, dotProduct_sum]
  exact Finset.sum_congr rfl fun i _ => by
    rw [Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul]

/-- The quadratic form of `secMat`. -/
theorem qform_secMat (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w : Fin M → ℝ) (γ : ℝ)
    (v : Fin r → ℝ) :
    v ⬝ᵥ (secMat θ R w γ *ᵥ v)
      = ∑ i, w i ^ 2 / (γ - w i ^ 2) * (v ⬝ᵥ (sigMat θ R i *ᵥ v)) := by
  simp only [secMat]
  exact qform_sum _ _ _

/-- The quadratic form of `secDerivMat`. -/
theorem qform_secDerivMat (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w : Fin M → ℝ) (γ : ℝ)
    (v : Fin r → ℝ) :
    v ⬝ᵥ (secDerivMat θ R w γ *ᵥ v)
      = ∑ i, w i ^ 2 / (γ - w i ^ 2) ^ 2 * (v ⬝ᵥ (sigMat θ R i *ᵥ v)) := by
  simp only [secDerivMat]
  exact qform_sum _ _ _

/-- The quadratic form of `swFmat2`. -/
theorem qform_swFmat2 (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (c w : Fin M → ℝ) (z : ℝ)
    (v : Fin r → ℝ) :
    v ⬝ᵥ (swFmat2 θ R c w z *ᵥ v)
      = (∑ i, w i ^ 2 * MPhet.ghetDeriv c w i z * (v ⬝ᵥ (sigMat θ R i *ᵥ v)))
        + MPhet.PsihetDeriv c w z * (v ⬝ᵥ v) := by
  simp only [swFmat2]
  rw [Matrix.add_mulVec, dotProduct_add, Matrix.smul_mulVec, Matrix.one_mulVec,
    dotProduct_smul, smul_eq_mul, qform_sum]

/-- `S_i = R_i Θ_i² R_iᵀ` is positive semidefinite. Private copy of the same fact in
`RankR/SingleWeight/Scalars.lean`, which is `private` there. -/
private theorem posSemidef_sigMat' (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (i : Fin M) :
    (sigMat θ R i).PosSemidef := by
  unfold sigMat
  rw [← Matrix.conjTranspose_eq_transpose_of_trivial (R i)]
  exact (Matrix.posSemidef_diagonal_iff.2 (fun j => sq_nonneg (θ i j))).mul_mul_conjTranspose_same
    (R i)

/-- Every `zᵀ S_i z` is nonnegative. -/
theorem qform_sigMat_nonneg (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (i : Fin M) (v : Fin r → ℝ) :
    0 ≤ v ⬝ᵥ (sigMat θ R i *ᵥ v) := by
  have h := (posSemidef_sigMat' θ R i).dotProduct_mulVec_nonneg v
  have hstar : star v = v := rfl
  rwa [hstar] at h

/-! ### 3. `phi` at `-1/γ`, and the outlier `swRho` -/

/-- `phi(-1/γ) = η(γ) = 1 - ∑_i c_i w_i⁴/(γ - w_i²)²` for every `γ` above `max_i w_i²`.
A copy of `MPhet.phi_neg_inv_gammaTop` (`RMT/Het/MPhet.lean:620`) with `gammaTop` replaced by
a bare `γ`. -/
theorem phi_neg_inv_eq {c w : Fin M → ℝ} {γ : ℝ} (hγ : Scalars.wSqMax w < γ) :
    MPhet.phi c w (-1 / γ) = 1 - ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 := by
  have hne : γ ≠ 0 := (gam_pos hγ).ne'
  unfold MPhet.phi
  congr 1
  refine Finset.sum_congr rfl fun i _ => ?_
  have hi := (gam_sub_pos hγ i).ne'
  have hden : 1 + w i ^ 2 * (-1 / γ) = (γ - w i ^ 2) / γ := by field_simp; ring
  rw [hden, div_pow, div_pow, div_div_eq_mul_div]
  field_simp

/-- `zfun'(-1/γ) = γ² η(γ)`, the twin of `MPhet.zfunDeriv_neg_inv_gammaTop`. -/
theorem zfunDeriv_neg_inv_eq {c w : Fin M → ℝ} {γ : ℝ} (hγ : Scalars.wSqMax w < γ) :
    MPhet.zfunDeriv c w (-1 / γ)
      = γ ^ 2 * (1 - ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2) := by
  have hγ0 : 0 < γ := gam_pos hγ
  have hne : γ ≠ 0 := hγ0.ne'
  have hne' : -1 / γ ≠ 0 := (div_neg_of_neg_of_pos (by norm_num) hγ0).ne
  have h1 := MPhet.phi_eq_mul c w hne'
  rw [phi_neg_inv_eq hγ] at h1
  rw [h1]
  field_simp

/-- Under the threshold `hη` the point `-1/γ` lies on the physical branch `(s⋆, 0)`. This is
the rank-`r` twin of the forward half of `MPhet.assumption4_iff_branch`. -/
theorem neg_inv_mem_Ioo {c w : Fin M → ℝ} {γ : ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hγ : Scalars.wSqMax w < γ) (hη : ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 < 1) :
    -1 / γ ∈ Set.Ioo (MPhet.sStar c w) 0 := by
  have hγ0 : 0 < γ := gam_pos hγ
  have hneg : -1 / γ < 0 := div_neg_of_neg_of_pos (by norm_num) hγ0
  have hW := MPhet.wSqMax_pos hw
  have hlo : MPhet.sLo w < -1 / γ := by
    have h := one_div_lt_one_div_of_lt hW hγ
    unfold MPhet.sLo
    rw [neg_div, neg_div]
    linarith
  have hmem : -1 / γ ∈ Set.Ioo (MPhet.sLo w) 0 := ⟨hlo, hneg⟩
  refine ⟨?_, hneg⟩
  have hphi : 0 < MPhet.phi c w (-1 / γ) := by
    rw [phi_neg_inv_eq hγ]; linarith
  by_contra hcon
  have hle : -1 / γ ≤ MPhet.sStar c w := not_lt.mp hcon
  have h := ((MPhet.phi_strictMonoOn hc hw).le_iff_le hmem (MPhet.sStar_mem hc hw)).mpr hle
  rw [MPhet.phi_sStar hc hw] at h
  linarith

/-- The outlier lies above the bulk edge. Rank-`r` twin of `MPhet.bHet_lt_rhoHet`
(`RMT/Het/MPhet.lean:673`). -/
theorem bHet_lt_swRho {c w : Fin M → ℝ} {γ : ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hγ : Scalars.wSqMax w < γ) (hη : ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 < 1) :
    MPhet.bHet c w < swRho c w γ := by
  have hmem := neg_inv_mem_Ioo hc hw hγ hη
  exact MPhet.bHet_lt_zfun hc hw hmem.1 hmem.2

/-- The outlier is positive: the bulk edge already is. -/
theorem swRho_pos {c w : Fin M → ℝ} {γ : ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hγ : Scalars.wSqMax w < γ) (hη : ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 < 1) :
    0 < swRho c w γ :=
  lt_trans (MPhet.bHet_pos hc hw) (bHet_lt_swRho hc hw hγ hη)

/-- The physical branch takes the outlier back to `-1/γ`. -/
theorem sPhys_swRho {c w : Fin M → ℝ} {γ : ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hγ : Scalars.wSqMax w < γ) (hη : ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 < 1) :
    MPhet.sPhys c w (swRho c w γ) = -1 / γ :=
  MPhet.sPhys_eq hc hw ⟨neg_inv_mem_Ioo hc hw hγ hη, rfl⟩

/-- `γ(ρ) = γ`: the secular variable of the outlier is the root it came from. Rank-`r` twin
of `MPhet.gamHet_rhoHet`. -/
theorem gamHet_swRho {c w : Fin M → ℝ} {γ : ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hγ : Scalars.wSqMax w < γ) (hη : ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 < 1) :
    MPhet.gamHet c w (swRho c w γ) = γ := by
  unfold MPhet.gamHet
  rw [sPhys_swRho hc hw hγ hη]
  have := (gam_pos hγ).ne'
  field_simp

/-- `s'(ρ) = 1/(γ² η(γ))`. Rank-`r` twin of `MPhet.sPhysDeriv_rhoHet`. -/
theorem sPhysDeriv_swRho {c w : Fin M → ℝ} {γ : ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hγ : Scalars.wSqMax w < γ) (hη : ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 < 1) :
    MPhet.sPhysDeriv c w (swRho c w γ)
      = 1 / (γ ^ 2 * (1 - ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2)) := by
  unfold MPhet.sPhysDeriv
  rw [sPhys_swRho hc hw hγ hη, zfunDeriv_neg_inv_eq hγ, one_div]

/-- `ρ = γ (1 + ∑_i c_i w_i²/(γ - w_i²))`, the closed form of the outlier. Rank-`r` twin of
`MPhet.rhoHet_eq_zfun` read backwards. -/
theorem swRho_eq_mul {c w : Fin M → ℝ} {γ : ℝ} (hγ : Scalars.wSqMax w < γ) :
    swRho c w γ = γ * (1 + ∑ i, c i * w i ^ 2 / (γ - w i ^ 2)) := by
  have hne : γ ≠ 0 := (gam_pos hγ).ne'
  unfold swRho MPhet.zfun
  rw [mul_add, mul_one, Finset.mul_sum]
  congr 1
  · field_simp
  · refine Finset.sum_congr rfl fun i _ => ?_
    have hi := (gam_sub_pos hγ i).ne'
    have hden : 1 + w i ^ 2 * (-1 / γ) = (γ - w i ^ 2) / γ := by field_simp; ring
    rw [hden, div_div_eq_mul_div]
    ring

/-- The outlier map is strictly increasing in the root. -/
theorem swRho_lt_swRho {c w : Fin M → ℝ} {γ γ' : ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hγ' : Scalars.wSqMax w < γ') (hη' : ∑ i, c i * w i ^ 4 / (γ' - w i ^ 2) ^ 2 < 1)
    (hγ : Scalars.wSqMax w < γ) (hη : ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 < 1)
    (h : γ' < γ) : swRho c w γ' < swRho c w γ := by
  have m' := neg_inv_mem_Ioo hc hw hγ' hη'
  have m := neg_inv_mem_Ioo hc hw hγ hη
  have hlt : -1 / γ' < -1 / γ := by
    rw [neg_div, neg_div, neg_lt_neg_iff]
    exact one_div_lt_one_div_of_lt (gam_pos hγ') h
  exact MPhet.zfun_strictMonoOn hc hw ⟨m'.1.le, m'.2⟩ ⟨m.1.le, m.2⟩ hlt

/-! ### 4. The matrix secular identity `I + F(z) = (γ/z)(I - secMat γ)` -/

/-- The scalar core of `one_add_swFmat_eq`: one entry of the matrix identity, with the
entry `(S_i)_{kl}` of the signal matrix and the entry `δ_{kl}` of the identity matrix left
as free reals. Mirrors the algebra of `MPhet.one_add_F_eq` (`RMT/Het/MPhet.lean:712`). -/
private theorem swFmat_scalar_aux {c w : Fin M → ℝ} {g : Fin M → ℝ} {z s : ℝ}
    (hz0 : 0 < z) (hs0 : s < 0) (hp : ∀ i, 0 < 1 + w i ^ 2 * s)
    (hzf : z = -1 / s + ∑ i, c i * w i ^ 2 / (1 + w i ^ 2 * s))
    (hg : ∀ i, g i = -1 / (z * (1 + w i ^ 2 * s)))
    (δ : ℝ) (S : Fin M → ℝ) :
    δ + ((∑ i, w i ^ 2 * g i * S i) + (∑ i, c i * w i ^ 2 * g i) * δ)
      = -1 / s / z * (δ - ∑ i, w i ^ 2 / (-1 / s - w i ^ 2) * S i) := by
  have hzne : z ≠ 0 := hz0.ne'
  have hsne : s ≠ 0 := hs0.ne
  have hpne : ∀ i, (1 + w i ^ 2 * s) ≠ 0 := fun i => (hp i).ne'
  have hden : ∀ i, -1 / s - w i ^ 2 = -(1 + w i ^ 2 * s) / s := by
    intro i; field_simp; ring
  have e1 : (∑ i, w i ^ 2 * g i * S i)
      = -(1 / z) * ∑ i, w i ^ 2 * S i / (1 + w i ^ 2 * s) := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    have hi := hpne i
    rw [hg i]
    field_simp
  have e2 : (∑ i, c i * w i ^ 2 * g i)
      = -(1 / z) * ∑ i, c i * w i ^ 2 / (1 + w i ^ 2 * s) := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    have hi := hpne i
    rw [hg i]
    field_simp
  have e3 : (∑ i, w i ^ 2 / (-1 / s - w i ^ 2) * S i)
      = -s * ∑ i, w i ^ 2 * S i / (1 + w i ^ 2 * s) := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    have hi := hpne i
    rw [hden i]
    field_simp
  have hB : (∑ i, c i * w i ^ 2 / (1 + w i ^ 2 * s)) = z + 1 / s := by
    have hneg : -1 / s = -(1 / s) := by ring
    rw [hneg] at hzf
    linarith
  rw [e1, e2, e3, hB]
  field_simp
  ring_nf

/-- **The matrix secular identity.** For every `z` above the bulk edge,
`I + F(z) = (γ(z)/z) (I - secMat θ R w γ(z))` with `γ(z) = -1/s(z)`. At `r = 1` this is
`MPhet.one_add_F_eq` (`RMT/Het/MPhet.lean:712`), where `det (I - secMat)` is
`Scalars.secular θ w γ`. Numeric check 1 of `notes/archive/trackG_specs/check_trackG.py`: max abs
error 3.3e-16 over 5 draws of `z`. -/
theorem one_add_swFmat_eq {θ : (i : Fin M) → Fin (rk i) → ℝ}
    {R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ} {c w : Fin M → ℝ}
    (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ} (hzb : MPhet.bHet c w < z) :
    1 + swFmat θ R c w z
      = (MPhet.gamHet c w z / z) • (1 - secMat θ R w (MPhet.gamHet c w z)) := by
  have hz0 : 0 < z := lt_trans (MPhet.bHet_pos hc hw) hzb
  have hs0 : MPhet.sPhys c w z < 0 := MPhet.sPhys_neg hc hw hzb
  have hp := MPhet.one_add_mul_sPhys_pos hc hw hzb
  have hzf : z = -1 / MPhet.sPhys c w z
      + ∑ i, c i * w i ^ 2 / (1 + w i ^ 2 * MPhet.sPhys c w z) := by
    have h := MPhet.zfun_sPhys hc hw hzb
    simp only [MPhet.zfun] at h
    exact h.symm
  have hg : ∀ i, MPhet.ghet c w i z
      = -1 / (z * (1 + w i ^ 2 * MPhet.sPhys c w z)) := fun _ => rfl
  ext k l
  simp only [swFmat, secMat, MPhet.Psihet, MPhet.gamHet, Matrix.add_apply, Matrix.sum_apply,
    Matrix.smul_apply, Matrix.sub_apply, smul_eq_mul]
  exact swFmat_scalar_aux hz0 hs0 hp hzf hg ((1 : Matrix (Fin r) (Fin r) ℝ) k l)
    (fun i => sigMat θ R i k l)

/-- **The outlier equation.** At `z = ρ` the eigenvector `z_ℓ` of `secMat` at eigenvalue `1`
satisfies `F(ρ) z_ℓ = -z_ℓ`, that is `(I + F(ρ)) z_ℓ = 0`. Numeric check 2 of
`notes/archive/trackG_specs/check_trackG.py`: max abs 1.3e-16 and 3.1e-16. -/
theorem swFmat_mulVec_eigvec {θ : (i : Fin M) → Fin (rk i) → ℝ}
    {R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ} {c w : Fin M → ℝ} {γ : ℝ}
    {z : EuclideanSpace ℝ (Fin r)} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hγ : Scalars.wSqMax w < γ) (hη : ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 < 1)
    (hz : secMat θ R w γ *ᵥ WithLp.ofLp z = WithLp.ofLp z) :
    swFmat θ R c w (swRho c w γ) *ᵥ WithLp.ofLp z = -WithLp.ofLp z := by
  have hid := one_add_swFmat_eq (θ := θ) (R := R) hc hw (bHet_lt_swRho hc hw hγ hη)
  rw [gamHet_swRho hc hw hγ hη] at hid
  have h0 := congrArg (fun A => A *ᵥ WithLp.ofLp z) hid
  simp only [Matrix.add_mulVec, Matrix.one_mulVec, Matrix.smul_mulVec, Matrix.sub_mulVec, hz,
    sub_self, smul_zero] at h0
  funext j
  have hj := congrFun h0 j
  simp only [Pi.add_apply, Pi.zero_apply] at hj
  simp only [Pi.neg_apply]
  linarith

/-! ### 5. The `ν` identity -/

/-- `g_i'(ρ)` in closed form: `γ/(ρ²(γ - w_i²)) + w_i²/(η(γ) ρ (γ - w_i²)²)`. This is the
termwise input of `qform_swFmat2_eq_swNu`, the twin of the `hT` step inside
`MPhet.rhoHet_mul_FhetDeriv` (`RMT/Het/MPhet.lean:788`). -/
theorem ghetDeriv_swRho {c w : Fin M → ℝ} {γ : ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hγ : Scalars.wSqMax w < γ) (hη : ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 < 1) (i : Fin M) :
    MPhet.ghetDeriv c w i (swRho c w γ)
      = γ / (swRho c w γ ^ 2 * (γ - w i ^ 2))
        + w i ^ 2 / ((1 - ∑ j, c j * w j ^ 4 / (γ - w j ^ 2) ^ 2) * swRho c w γ
            * (γ - w i ^ 2) ^ 2) := by
  have hγ0 : 0 < γ := gam_pos hγ
  have h1 : (γ - w i ^ 2) ≠ 0 := (gam_sub_pos hγ i).ne'
  have h2 : γ ≠ 0 := hγ0.ne'
  have h3 : (1 - ∑ j, c j * w j ^ 4 / (γ - w j ^ 2) ^ 2) ≠ 0 := by
    have hpos : 0 < 1 - ∑ j, c j * w j ^ 4 / (γ - w j ^ 2) ^ 2 := by linarith
    exact hpos.ne'
  have h4 : swRho c w γ ≠ 0 := (swRho_pos hc hw hγ hη).ne'
  simp only [MPhet.ghetDeriv]
  rw [sPhys_swRho hc hw hγ hη, sPhysDeriv_swRho hc hw hγ hη]
  have hden : 1 + w i ^ 2 * (-1 / γ) = (γ - w i ^ 2) / γ := by field_simp; ring
  rw [hden]
  field_simp

/-- **The `ν` identity.** `z_ℓᵀ F'(ρ_ℓ) z_ℓ = γ_ℓ (z_ℓᵀ K(γ_ℓ) z_ℓ)/(η(γ_ℓ) ρ_ℓ)`, the
rank-`r` twin of `MPhet.rhoHet_mul_FhetDeriv` (`RMT/Het/MPhet.lean:788`). Numeric check 3b of
`notes/archive/trackG_specs/check_trackG.py`: `|diff|` 1.6e-10 and 2.0e-10, the residual being the
central-difference truncation of the numeric `F'`. -/
theorem qform_swFmat2_eq_swNu {θ : (i : Fin M) → Fin (rk i) → ℝ}
    {R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ} {c w : Fin M → ℝ} {γ : ℝ}
    {z : EuclideanSpace ℝ (Fin r)} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hγ : Scalars.wSqMax w < γ) (hη : ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 < 1)
    (hz : secMat θ R w γ *ᵥ WithLp.ofLp z = WithLp.ofLp z) :
    WithLp.ofLp z ⬝ᵥ (swFmat2 θ R c w (swRho c w γ) *ᵥ WithLp.ofLp z)
      = swNu θ R c w γ z := by
  have hγ0 : 0 < γ := gam_pos hγ
  have hη0 : 0 < 1 - ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 := by linarith
  have hρ0 : 0 < swRho c w γ := swRho_pos hc hw hγ hη
  have h2 : γ ≠ 0 := hγ0.ne'
  have h3 : (1 - ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2) ≠ 0 := hη0.ne'
  have h4 : swRho c w γ ≠ 0 := hρ0.ne'
  have hvv : WithLp.ofLp z ⬝ᵥ WithLp.ofLp z
      = ∑ i, w i ^ 2 / (γ - w i ^ 2)
          * (WithLp.ofLp z ⬝ᵥ (sigMat θ R i *ᵥ WithLp.ofLp z)) := by
    have h := qform_secMat θ R w γ (WithLp.ofLp z)
    rw [hz] at h
    exact h
  have hA : (∑ i, w i ^ 2 * MPhet.ghetDeriv c w i (swRho c w γ)
        * (WithLp.ofLp z ⬝ᵥ (sigMat θ R i *ᵥ WithLp.ofLp z)))
      = γ / swRho c w γ ^ 2
          * (∑ i, w i ^ 2 / (γ - w i ^ 2)
              * (WithLp.ofLp z ⬝ᵥ (sigMat θ R i *ᵥ WithLp.ofLp z)))
        + 1 / ((1 - ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2) * swRho c w γ)
          * (∑ i, w i ^ 4 / (γ - w i ^ 2) ^ 2
              * (WithLp.ofLp z ⬝ᵥ (sigMat θ R i *ᵥ WithLp.ofLp z))) := by
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    have h1 : (γ - w i ^ 2) ≠ 0 := (gam_sub_pos hγ i).ne'
    rw [ghetDeriv_swRho hc hw hγ hη i]
    field_simp
  have hB : MPhet.PsihetDeriv c w (swRho c w γ)
      = γ / swRho c w γ ^ 2 * (∑ i, c i * w i ^ 2 / (γ - w i ^ 2))
        + 1 / ((1 - ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2) * swRho c w γ)
          * (∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2) := by
    simp only [MPhet.PsihetDeriv]
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    have h1 : (γ - w i ^ 2) ≠ 0 := (gam_sub_pos hγ i).ne'
    rw [ghetDeriv_swRho hc hw hγ hη i]
    field_simp
  have hC : (∑ i, c i * w i ^ 2 / (γ - w i ^ 2)) = swRho c w γ / γ - 1 := by
    rw [eq_sub_iff_add_eq, eq_div_iff h2, swRho_eq_mul hγ]
    ring
  have hDN : (∑ i, w i ^ 4 / (γ - w i ^ 2) ^ 2
        * (WithLp.ofLp z ⬝ᵥ (sigMat θ R i *ᵥ WithLp.ofLp z)))
      + (∑ i, w i ^ 2 / (γ - w i ^ 2)
        * (WithLp.ofLp z ⬝ᵥ (sigMat θ R i *ᵥ WithLp.ofLp z)))
      = γ * ∑ i, w i ^ 2 / (γ - w i ^ 2) ^ 2
          * (WithLp.ofLp z ⬝ᵥ (sigMat θ R i *ᵥ WithLp.ofLp z)) := by
    rw [← Finset.sum_add_distrib, Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    have h1 : (γ - w i ^ 2) ≠ 0 := (gam_sub_pos hγ i).ne'
    field_simp
    ring
  rw [qform_swFmat2, hA, hB, hvv, hC]
  simp only [swNu]
  rw [qform_secDerivMat, ← hDN]
  field_simp
  ring

/-- `zᵀ K(γ) z > 0` for a unit eigenvector: termwise `w_i²/(γ - w_i²)² ≥ (1/γ) w_i²/(γ - w_i²)`
and `zᵀ secMat(γ) z = ‖z‖² = 1`, so the form is at least `1/γ`. -/
theorem qform_secDerivMat_pos {θ : (i : Fin M) → Fin (rk i) → ℝ}
    {R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ} {w : Fin M → ℝ} {γ : ℝ}
    {z : EuclideanSpace ℝ (Fin r)} (hγ : Scalars.wSqMax w < γ) (hznorm : ‖z‖ = 1)
    (hz : secMat θ R w γ *ᵥ WithLp.ofLp z = WithLp.ofLp z) :
    0 < WithLp.ofLp z ⬝ᵥ (secDerivMat θ R w γ *ᵥ WithLp.ofLp z) := by
  have hγ0 : 0 < γ := gam_pos hγ
  have hvv : WithLp.ofLp z ⬝ᵥ WithLp.ofLp z = 1 := by
    have h := real_inner_eq_dotProduct z z
    rw [real_inner_self_eq_norm_sq, hznorm] at h
    simpa using h.symm
  have h1 : (∑ i, w i ^ 2 / (γ - w i ^ 2)
      * (WithLp.ofLp z ⬝ᵥ (sigMat θ R i *ᵥ WithLp.ofLp z))) = 1 := by
    have h := qform_secMat θ R w γ (WithLp.ofLp z)
    rw [hz, hvv] at h
    exact h.symm
  have h2 : 1 / γ * (∑ i, w i ^ 2 / (γ - w i ^ 2)
        * (WithLp.ofLp z ⬝ᵥ (sigMat θ R i *ᵥ WithLp.ofLp z)))
      ≤ ∑ i, w i ^ 2 / (γ - w i ^ 2) ^ 2
        * (WithLp.ofLp z ⬝ᵥ (sigMat θ R i *ᵥ WithLp.ofLp z)) := by
    rw [Finset.mul_sum]
    refine Finset.sum_le_sum fun i _ => ?_
    have hQ := qform_sigMat_nonneg θ R i (WithLp.ofLp z)
    have hD : 0 < γ - w i ^ 2 := gam_sub_pos hγ i
    have hDne : (γ - w i ^ 2) ≠ 0 := hD.ne'
    have hγne : γ ≠ 0 := hγ0.ne'
    have e : w i ^ 2 / (γ - w i ^ 2) ^ 2 - 1 / γ * (w i ^ 2 / (γ - w i ^ 2))
        = w i ^ 4 / (γ * (γ - w i ^ 2) ^ 2) := by
      field_simp
      ring
    have hnn : 0 ≤ w i ^ 4 / (γ * (γ - w i ^ 2) ^ 2) := by positivity
    have hcoef : 1 / γ * (w i ^ 2 / (γ - w i ^ 2)) ≤ w i ^ 2 / (γ - w i ^ 2) ^ 2 := by
      linarith
    calc 1 / γ * (w i ^ 2 / (γ - w i ^ 2)
            * (WithLp.ofLp z ⬝ᵥ (sigMat θ R i *ᵥ WithLp.ofLp z)))
        = 1 / γ * (w i ^ 2 / (γ - w i ^ 2))
            * (WithLp.ofLp z ⬝ᵥ (sigMat θ R i *ᵥ WithLp.ofLp z)) := by ring
      _ ≤ w i ^ 2 / (γ - w i ^ 2) ^ 2
            * (WithLp.ofLp z ⬝ᵥ (sigMat θ R i *ᵥ WithLp.ofLp z)) :=
          mul_le_mul_of_nonneg_right hcoef hQ
  rw [h1] at h2
  have hpos : 0 < 1 / γ := by positivity
  rw [qform_secDerivMat]
  linarith

/-- `ν_ℓ > 0`. -/
theorem swNu_pos {θ : (i : Fin M) → Fin (rk i) → ℝ}
    {R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ} {c w : Fin M → ℝ} {γ : ℝ}
    {z : EuclideanSpace ℝ (Fin r)} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hγ : Scalars.wSqMax w < γ) (hη : ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 < 1)
    (hznorm : ‖z‖ = 1) (hz : secMat θ R w γ *ᵥ WithLp.ofLp z = WithLp.ofLp z) :
    0 < swNu θ R c w γ z := by
  have hγ0 : 0 < γ := gam_pos hγ
  have hη0 : 0 < 1 - ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 := by linarith
  have hρ0 : 0 < swRho c w γ := swRho_pos hc hw hγ hη
  have hK := qform_secDerivMat_pos hγ hznorm hz
  unfold swNu
  exact div_pos (mul_pos hγ0 hK) (mul_pos hη0 hρ0)

/-- **The overlap identity.** `1/(ρ_ℓ ν_ℓ) = swTerm θ R w c γ_ℓ z_ℓ`, the summand of
`main_paper.tex:2119`. Rank-`r` twin of `MPhet.overlap_identity_het`
(`RMT/Het/MPhet.lean:842`). Numeric check 3 of `notes/archive/trackG_specs/check_trackG.py`:
0.8869271 versus 0.8869271 and 0.8505194 versus 0.8505194. -/
theorem one_div_swRho_mul_swNu {θ : (i : Fin M) → Fin (rk i) → ℝ}
    {R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ} {c w : Fin M → ℝ} {γ : ℝ}
    {z : EuclideanSpace ℝ (Fin r)} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hγ : Scalars.wSqMax w < γ) (hη : ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 < 1)
    (hznorm : ‖z‖ = 1) (hz : secMat θ R w γ *ᵥ WithLp.ofLp z = WithLp.ofLp z) :
    1 / (swRho c w γ * swNu θ R c w γ z) = swTerm θ R w c γ z := by
  have hγ0 : 0 < γ := gam_pos hγ
  have hη0 : 0 < 1 - ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 := by linarith
  have hρ0 : 0 < swRho c w γ := swRho_pos hc hw hγ hη
  have hK := qform_secDerivMat_pos hγ hznorm hz
  have h1 : γ ≠ 0 := hγ0.ne'
  have h2 : (1 - ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2) ≠ 0 := hη0.ne'
  have h3 : swRho c w γ ≠ 0 := hρ0.ne'
  have h4 : (WithLp.ofLp z ⬝ᵥ (secDerivMat θ R w γ *ᵥ WithLp.ofLp z)) ≠ 0 := hK.ne'
  simp only [swNu, swTerm]
  field_simp

end SingleWeight

end StackedSVD

-- Deciding check of unit G0 (`notes/archive/trackG_plan.md`), both elaborated on 2026-09-05:
-- #check @StackedSVD.SingleWeight.one_add_swFmat_eq
-- #check @StackedSVD.SingleWeight.one_div_swRho_mul_swNu
