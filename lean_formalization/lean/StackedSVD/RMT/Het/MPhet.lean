/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.StackSVDWeighted
import StackedSVD.RMT.MP

/-!
# Item MP-het: the scalar layer of the heteroscedastic Marchenko-Pastur law

Task H1 of `notes/archive/plan_heterolaw_A.md` (sections 2.2 to 2.4 and 3.5). The weighted stack
`X_stack(w) = [w_1 X_1; ...; w_M X_M]` of `main_paper.tex:1350` has block heteroscedastic
noise `Σ^{1/2} E`, `Σ = diag(w_i² I_{n_i})`. Its noise Gram matrix has the Silverstein law:
the `d`-side trace `s(z)` solves `z = zfun s := -1/s + ∑_i c_i w_i²/(1 + w_i² s)`. Everything
here is real algebra on the physical branch of this equation, in the style of `RMT/MP.lean`.

## Objects

* `sLo w = -1/max_i w_i²`, the left end of the domain of `zfun`.
* `phi c w s = s² zfun'(s) = 1 - ∑ c_i w_i⁴ s²/(1 + w_i² s)²`, strictly increasing on
  `(sLo, 0)` from `-∞` to `1`; `sStar c w` is its unique zero there.
* `bHet c w = zfun (sStar)`, the upper edge of the bulk. `zfun` decreases on `(sLo, sStar]`
  and increases on `[sStar, 0)` to `+∞`.
* `sPhys c w z`, the inverse of `zfun` on `(sStar, 0)`, defined for `z > bHet`; it is
  strictly increasing, differentiable with `sPhysDeriv = 1/zfun'(sPhys z)`, and
  `sPhysDeriv → +∞` at `bHet⁺` (the square-root edge, the R3⁻ input).
* `rhoHet θ c w = γ₁ (1 + ∑ c_i w_i²/(γ₁ - w_i²))`, `γ₁ = gammaTop θ w`: the outlier.
* `ghet c w i z = -1/(z (1 + w_i² sPhys z))`, the block-`i` partial trace; `Phihet`,
  `Psihet`, `Fhet = Phihet + Psihet` the limiting secular function `1 + F(z)` of section 2.4.
* `bSF c w = (max_i |w_i| + √(∑ c_i w_i²))²`, the Sudakov-Fernique edge of section 3.6.

## Results

`assumption4_iff_branch` (`eq:assumption4`, `main_paper.tex:1368`, is exactly the branch
condition `-1/γ₁ ∈ (sStar, 0)`), `bHet_lt_rhoHet`, `one_add_F_eq` (`1 + F(z) = (γ/z)
secular θ w γ` with `γ = -1/sPhys z`), `one_add_F_rhoHet` and `one_add_F_eq_zero_iff` (the
outlier is the unique root above the edge), `overlap_identity_het` (`1/(ρ F'(ρ)) = L(w)`,
the paper's `L(w)` of `main_paper.tex:1433`), `bHet_le_bSF`, and the `M = 1` reduction to
`RMT/MP.lean` (`sStar_single`, `bHet_single`, `sPhys_single`, `rhoHet_single`).

## Junk conventions

`sLo 0 = 0`; `sStar` and `sPhys` are `0` when the defining root does not exist or is not
unique (never the case under `0 < c i` and some `w i ≠ 0`, respectively `bHet < z`);
`gammaTop = 0` when no root exists, so `rhoHet = 0` there. Every lemma carries the
hypotheses that exclude the junk.

Numeric check of every identity: a session script (`check_mphet.py`, seed 20260901; not
kept).
-/

open Filter Topology Finset

namespace StackedSVD
namespace MPhet

open Scalars

variable {M : ℕ}

/-! ### Definitions -/

/-- `-1/max_i w_i²`, the left end of the domain of `zfun`. Junk `-1/0 = 0` at `w = 0`. -/
noncomputable def sLo (w : Fin M → ℝ) : ℝ := -1 / wSqMax w

/-- The Silverstein map `zfun s = -1/s + ∑ c_i w_i²/(1 + w_i² s)` (plan section 2.2). -/
noncomputable def zfun (c w : Fin M → ℝ) (s : ℝ) : ℝ :=
  -1 / s + ∑ i, c i * w i ^ 2 / (1 + w i ^ 2 * s)

/-- `zfun'(s) = 1/s² - ∑ c_i w_i⁴/(1 + w_i² s)²`. -/
noncomputable def zfunDeriv (c w : Fin M → ℝ) (s : ℝ) : ℝ :=
  1 / s ^ 2 - ∑ i, c i * w i ^ 4 / (1 + w i ^ 2 * s) ^ 2

/-- `phi s = s² zfun'(s) = 1 - ∑ c_i w_i⁴ s²/(1 + w_i² s)²` (plan section 2.3). -/
noncomputable def phi (c w : Fin M → ℝ) (s : ℝ) : ℝ :=
  1 - ∑ i, c i * w i ^ 4 * s ^ 2 / (1 + w i ^ 2 * s) ^ 2

/-- `s` is the critical point: a zero of `phi` in `(sLo, 0)`. -/
def IsSStar (c w : Fin M → ℝ) (s : ℝ) : Prop := s ∈ Set.Ioo (sLo w) 0 ∧ phi c w s = 0

open Classical in
/-- The critical point `s⋆`, the unique zero of `phi` in `(sLo, 0)`. Junk `0` otherwise. -/
noncomputable def sStar (c w : Fin M → ℝ) : ℝ :=
  if h : ∃! s, IsSStar c w s then h.choose else 0

/-- The upper bulk edge `bHet = zfun (s⋆)` of the noise Gram matrix (plan section 2.3). -/
noncomputable def bHet (c w : Fin M → ℝ) : ℝ := zfun c w (sStar c w)

/-- `s` is the physical root: `zfun s = z` with `s ∈ (s⋆, 0)`. -/
def IsSPhys (c w : Fin M → ℝ) (z s : ℝ) : Prop := s ∈ Set.Ioo (sStar c w) 0 ∧ zfun c w s = z

open Classical in
/-- The physical branch `s(z)`, the inverse of `zfun` on `(s⋆, 0)`. Junk `0` for `z ≤ bHet`. -/
noncomputable def sPhys (c w : Fin M → ℝ) (z : ℝ) : ℝ :=
  if h : ∃! s, IsSPhys c w z s then h.choose else 0

/-- `s'(z) = 1/zfun'(s(z))`, the derivative of the physical branch. -/
noncomputable def sPhysDeriv (c w : Fin M → ℝ) (z : ℝ) : ℝ :=
  (zfunDeriv c w (sPhys c w z))⁻¹

/-- `γ(z) = -1/s(z)`, the variable in which `1 + F` becomes the paper's secular function. -/
noncomputable def gamHet (c w : Fin M → ℝ) (z : ℝ) : ℝ := -1 / sPhys c w z

/-- The outlier `ρ = γ₁ (1 + ∑ c_i w_i²/(γ₁ - w_i²)) = zfun (-1/γ₁)` (plan section 2.4).
Junk `0` when `gammaTop θ w = 0`. -/
noncomputable def rhoHet (θ c w : Fin M → ℝ) : ℝ :=
  gammaTop θ w * (1 + ∑ i, c i * w i ^ 2 / (gammaTop θ w - w i ^ 2))

/-- The Sudakov-Fernique edge `(max_i |w_i| + √(∑ c_i w_i²))²` (plan section 3.6), written
with `√(wSqMax w) = max_i |w_i|`. -/
noncomputable def bSF (c w : Fin M → ℝ) : ℝ :=
  (Real.sqrt (wSqMax w) + Real.sqrt (∑ i, c i * w i ^ 2)) ^ 2

/-- The block-`i` partial trace `g_i(z) = -1/(z (1 + w_i² s(z)))` (plan section 2.2). -/
noncomputable def ghet (c w : Fin M → ℝ) (i : Fin M) (z : ℝ) : ℝ :=
  -1 / (z * (1 + w i ^ 2 * sPhys c w z))

/-- `g_i'(z)`, in closed form through `sPhysDeriv`. -/
noncomputable def ghetDeriv (c w : Fin M → ℝ) (i : Fin M) (z : ℝ) : ℝ :=
  (1 + w i ^ 2 * sPhys c w z + z * w i ^ 2 * sPhysDeriv c w z)
    / (z * (1 + w i ^ 2 * sPhys c w z)) ^ 2

/-- `Φ(z) = ∑ θ_i² w_i² g_i(z)`, the limit of `ũ₀ᵀ G₀' ũ₀`. -/
noncomputable def Phihet (θ c w : Fin M → ℝ) (z : ℝ) : ℝ :=
  ∑ i, θ i ^ 2 * w i ^ 2 * ghet c w i z

/-- `Ψ(z) = ∑ c_i w_i² g_i(z)`, the limit of `eᵀ Σ^{1/2} G₀' Σ^{1/2} e`. -/
noncomputable def Psihet (c w : Fin M → ℝ) (z : ℝ) : ℝ := ∑ i, c i * w i ^ 2 * ghet c w i z

/-- `F = Φ + Ψ`; the limiting secular function is `1 + F`. -/
noncomputable def Fhet (θ c w : Fin M → ℝ) (z : ℝ) : ℝ := Phihet θ c w z + Psihet c w z

/-- `Φ'(z)`, in closed form. -/
noncomputable def PhihetDeriv (θ c w : Fin M → ℝ) (z : ℝ) : ℝ :=
  ∑ i, θ i ^ 2 * w i ^ 2 * ghetDeriv c w i z

/-- `Ψ'(z)`, in closed form. -/
noncomputable def PsihetDeriv (c w : Fin M → ℝ) (z : ℝ) : ℝ :=
  ∑ i, c i * w i ^ 2 * ghetDeriv c w i z

/-- `F'(z) = Φ'(z) + Ψ'(z)`. -/
noncomputable def FhetDeriv (θ c w : Fin M → ℝ) (z : ℝ) : ℝ :=
  PhihetDeriv θ c w z + PsihetDeriv c w z

/-! ### The domain `(sLo, 0)` -/

theorem wSqMax_pos {w : Fin M → ℝ} (hw : ∃ i, w i ≠ 0) : 0 < wSqMax w := by
  obtain ⟨i, hi⟩ := hw
  exact lt_of_lt_of_le (by positivity) (le_wSqMax w i)

theorem sLo_neg {w : Fin M → ℝ} (hw : ∃ i, w i ≠ 0) : sLo w < 0 :=
  div_neg_of_neg_of_pos (by norm_num) (wSqMax_pos hw)

/-- For `s > sLo`, every `1 + w_i² s` is positive. -/
theorem one_add_mul_pos {w : Fin M → ℝ} (hw : ∃ i, w i ≠ 0) {s : ℝ} (hs : sLo w < s)
    (i : Fin M) : 0 < 1 + w i ^ 2 * s := by
  have hW := wSqMax_pos hw
  have hle := le_wSqMax w i
  rcases le_or_gt 0 s with h0 | h0
  · have : 0 ≤ w i ^ 2 * s := mul_nonneg (sq_nonneg _) h0
    linarith
  · have h1 : -1 < s * wSqMax w := by
      unfold sLo at hs
      exact (div_lt_iff₀ hW).mp hs
    have h2 : wSqMax w * s ≤ w i ^ 2 * s := mul_le_mul_of_nonpos_right hle h0.le
    nlinarith

/-! ### `phi`: strictly increasing on `(sLo, 0)`, with a zero `sStar` -/

/-- Each summand of `phi` is antitone on the domain. -/
theorem phi_term_le {c w : Fin M → ℝ} (hw : ∃ i, w i ≠ 0) (i : Fin M) (hci : 0 ≤ c i)
    {s₁ s₂ : ℝ} (h₁ : sLo w < s₁) (h₂ : s₂ < 0) (hlt : s₁ < s₂) :
    c i * w i ^ 4 * s₂ ^ 2 / (1 + w i ^ 2 * s₂) ^ 2
      ≤ c i * w i ^ 4 * s₁ ^ 2 / (1 + w i ^ 2 * s₁) ^ 2 := by
  have p₁ := one_add_mul_pos hw h₁ i
  have p₂ := one_add_mul_pos hw (h₁.trans hlt) i
  rw [div_le_div_iff₀ (by positivity) (by positivity)]
  have hyx : w i ^ 2 * s₁ * (1 + w i ^ 2 * s₂) ≤ w i ^ 2 * s₂ * (1 + w i ^ 2 * s₁) := by
    nlinarith [mul_nonneg (sq_nonneg (w i)) (sub_nonneg.2 hlt.le)]
  have hws : w i ^ 2 * s₂ ≤ 0 := by nlinarith [sq_nonneg (w i)]
  have hx0 : w i ^ 2 * s₂ * (1 + w i ^ 2 * s₁) ≤ 0 := by nlinarith
  have hsq : (w i ^ 2 * s₂ * (1 + w i ^ 2 * s₁)) ^ 2
      ≤ (w i ^ 2 * s₁ * (1 + w i ^ 2 * s₂)) ^ 2 := by
    rw [← neg_sq (w i ^ 2 * s₁ * (1 + w i ^ 2 * s₂))]
    exact sq_le_sq' (by linarith) (by linarith)
  calc c i * w i ^ 4 * s₂ ^ 2 * (1 + w i ^ 2 * s₁) ^ 2
      = c i * (w i ^ 2 * s₂ * (1 + w i ^ 2 * s₁)) ^ 2 := by ring
    _ ≤ c i * (w i ^ 2 * s₁ * (1 + w i ^ 2 * s₂)) ^ 2 := mul_le_mul_of_nonneg_left hsq hci
    _ = c i * w i ^ 4 * s₁ ^ 2 * (1 + w i ^ 2 * s₂) ^ 2 := by ring

/-- The summand of a table with `w i ≠ 0` and `0 < c i` is strictly antitone. -/
theorem phi_term_lt {c w : Fin M → ℝ} (hw : ∃ i, w i ≠ 0) (i : Fin M) (hwi : w i ≠ 0)
    (hci : 0 < c i) {s₁ s₂ : ℝ} (h₁ : sLo w < s₁) (h₂ : s₂ < 0) (hlt : s₁ < s₂) :
    c i * w i ^ 4 * s₂ ^ 2 / (1 + w i ^ 2 * s₂) ^ 2
      < c i * w i ^ 4 * s₁ ^ 2 / (1 + w i ^ 2 * s₁) ^ 2 := by
  have p₁ := one_add_mul_pos hw h₁ i
  have p₂ := one_add_mul_pos hw (h₁.trans hlt) i
  have hw2 : 0 < w i ^ 2 := by positivity
  rw [div_lt_div_iff₀ (by positivity) (by positivity)]
  have hyx : w i ^ 2 * s₁ * (1 + w i ^ 2 * s₂) < w i ^ 2 * s₂ * (1 + w i ^ 2 * s₁) := by
    nlinarith [mul_pos hw2 (sub_pos.2 hlt)]
  have hws : w i ^ 2 * s₂ < 0 := mul_neg_of_pos_of_neg hw2 h₂
  have hx0 : w i ^ 2 * s₂ * (1 + w i ^ 2 * s₁) < 0 := mul_neg_of_neg_of_pos hws p₁
  have hsq : (w i ^ 2 * s₂ * (1 + w i ^ 2 * s₁)) ^ 2
      < (w i ^ 2 * s₁ * (1 + w i ^ 2 * s₂)) ^ 2 := by
    rw [← neg_sq (w i ^ 2 * s₁ * (1 + w i ^ 2 * s₂))]
    exact sq_lt_sq' (by linarith) (by linarith)
  calc c i * w i ^ 4 * s₂ ^ 2 * (1 + w i ^ 2 * s₁) ^ 2
      = c i * (w i ^ 2 * s₂ * (1 + w i ^ 2 * s₁)) ^ 2 := by ring
    _ < c i * (w i ^ 2 * s₁ * (1 + w i ^ 2 * s₂)) ^ 2 := mul_lt_mul_of_pos_left hsq hci
    _ = c i * w i ^ 4 * s₁ ^ 2 * (1 + w i ^ 2 * s₂) ^ 2 := by ring

theorem phi_strictMonoOn {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    StrictMonoOn (phi c w) (Set.Ioo (sLo w) 0) := by
  intro s₁ hs₁ s₂ hs₂ hlt
  obtain ⟨k, hk⟩ := hw
  have hsum : ∑ i, c i * w i ^ 4 * s₂ ^ 2 / (1 + w i ^ 2 * s₂) ^ 2
      < ∑ i, c i * w i ^ 4 * s₁ ^ 2 / (1 + w i ^ 2 * s₁) ^ 2 :=
    Finset.sum_lt_sum (fun i _ => phi_term_le ⟨k, hk⟩ i (hc i).le hs₁.1 hs₂.2 hlt)
      ⟨k, Finset.mem_univ _, phi_term_lt ⟨k, hk⟩ k hk (hc k) hs₁.1 hs₂.2 hlt⟩
  unfold phi
  linarith

theorem continuousOn_phi {c w : Fin M → ℝ} (hw : ∃ i, w i ≠ 0) :
    ContinuousOn (phi c w) (Set.Ioi (sLo w)) := by
  unfold phi
  refine continuousOn_const.sub (continuousOn_finsetSum _ fun i _ => ?_)
  refine ContinuousOn.div (by fun_prop) (by fun_prop) fun s hs => ?_
  exact pow_ne_zero 2 (one_add_mul_pos hw hs i).ne'

theorem phi_zero (c w : Fin M → ℝ) : phi c w 0 = 1 := by simp [phi]

/-- A point near `0` where `phi` is positive. -/
theorem exists_phi_pos {c w : Fin M → ℝ} (hw : ∃ i, w i ≠ 0) :
    ∃ b ∈ Set.Ioo (sLo w) 0, 0 < phi c w b := by
  have hcont : ContinuousAt (phi c w) 0 :=
    (continuousOn_phi hw).continuousAt (Ioi_mem_nhds (sLo_neg hw))
  have h1 : ∀ᶠ s in 𝓝 (0:ℝ), 0 < phi c w s := by
    have : (0:ℝ) < phi c w 0 := by rw [phi_zero]; norm_num
    exact hcont.eventually (lt_mem_nhds this)
  have h2 : ∀ᶠ s in 𝓝[<] (0:ℝ), 0 < phi c w s := h1.filter_mono nhdsWithin_le_nhds
  have h3 : ∀ᶠ s in 𝓝[<] (0:ℝ), s ∈ Set.Ioo (sLo w) 0 := Ioo_mem_nhdsLT (sLo_neg hw)
  obtain ⟨b, hb1, hb2⟩ := (h3.and h2).exists
  exact ⟨b, hb1, hb2⟩

/-- A point near `sLo` where `phi` is negative: the summand of a table `k` that attains
`max_i w_i²` alone exceeds `1` at `s = -(1 - δ)/max_i w_i²` with `δ = √c_k/(2(1 + √c_k))`. -/
theorem exists_phi_neg {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    ∃ a ∈ Set.Ioo (sLo w) 0, phi c w a < 0 := by
  have : Nonempty (Fin M) := ⟨hw.choose⟩
  obtain ⟨k, hk⟩ := exists_eq_ciSup_of_finite (f := fun i => w i ^ 2)
  have hk' : w k ^ 2 = wSqMax w := hk
  have hW := wSqMax_pos hw
  set r := Real.sqrt (c k) with hr
  have hr0 : 0 < r := Real.sqrt_pos.mpr (hc k)
  have hr2 : r ^ 2 = c k := Real.sq_sqrt (hc k).le
  set δ := r / (2 * (1 + r)) with hδ
  have hδ0 : 0 < δ := by positivity
  have hδ1 : δ < 1 / 2 := by
    rw [hδ, div_lt_div_iff₀ (by positivity) (by positivity)]; nlinarith
  have hδr : δ < r / 2 := by
    rw [hδ, div_lt_div_iff₀ (by positivity) (by positivity)]; nlinarith
  refine ⟨-(1 - δ) / wSqMax w, ⟨?_, ?_⟩, ?_⟩
  · unfold sLo; rw [div_lt_div_iff_of_pos_right hW]; linarith
  · exact div_neg_of_neg_of_pos (by linarith) hW
  · unfold phi
    have hterm : ∀ i, 0 ≤ c i * w i ^ 4 * (-(1 - δ) / wSqMax w) ^ 2
        / (1 + w i ^ 2 * (-(1 - δ) / wSqMax w)) ^ 2 := fun i => by
      have := (hc i).le; positivity
    have hk_term : c k * w k ^ 4 * (-(1 - δ) / wSqMax w) ^ 2
        / (1 + w k ^ 2 * (-(1 - δ) / wSqMax w)) ^ 2 = c k * (1 - δ) ^ 2 / δ ^ 2 := by
      rw [show w k ^ 4 = (w k ^ 2) ^ 2 by ring, hk']
      have hδne : δ ≠ 0 := hδ0.ne'
      have hden : 1 + wSqMax w * (-(1 - δ) / wSqMax w) = δ := by field_simp; ring
      rw [hden]
      field_simp
    have hle : c k * (1 - δ) ^ 2 / δ ^ 2 ≤ ∑ i, c i * w i ^ 4 * (-(1 - δ) / wSqMax w) ^ 2
        / (1 + w i ^ 2 * (-(1 - δ) / wSqMax w)) ^ 2 := by
      rw [← hk_term]
      exact Finset.single_le_sum (fun i _ => hterm i) (Finset.mem_univ k)
    have hgt : 1 < c k * (1 - δ) ^ 2 / δ ^ 2 := by
      rw [lt_div_iff₀ (by positivity)]
      have h1 : δ ^ 2 < c k / 4 := by
        nlinarith [mul_pos (sub_pos.2 hδr) (add_pos (half_pos hr0) hδ0)]
      have h2 : 1 / 4 < (1 - δ) ^ 2 := by
        nlinarith [mul_pos (sub_pos.2 hδ1) (show (0:ℝ) < 3 / 2 - δ by linarith)]
      have h3 : c k / 4 < c k * (1 - δ) ^ 2 := by
        have := mul_lt_mul_of_pos_left h2 (hc k); linarith
      linarith
    linarith

theorem exists_isSStar {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    ∃ s, IsSStar c w s := by
  obtain ⟨a, ha, hfa⟩ := exists_phi_neg hc hw
  obtain ⟨b, hb, hfb⟩ := exists_phi_pos (c := c) hw
  have hab : a < b := ((phi_strictMonoOn hc hw).lt_iff_lt ha hb).mp (by linarith)
  have hcont : ContinuousOn (phi c w) (Set.Icc a b) :=
    (continuousOn_phi hw).mono fun s hs => lt_of_lt_of_le ha.1 hs.1
  obtain ⟨s, hs, hfs⟩ := intermediate_value_Ioo hab.le hcont ⟨hfa, hfb⟩
  exact ⟨s, ⟨lt_trans ha.1 hs.1, lt_trans hs.2 hb.2⟩, hfs⟩

theorem isSStar_unique {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {s t : ℝ}
    (hs : IsSStar c w s) (ht : IsSStar c w t) : s = t :=
  (phi_strictMonoOn hc hw).injOn hs.1 ht.1 (hs.2.trans ht.2.symm)

theorem existsUnique_isSStar {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    ∃! s, IsSStar c w s := by
  obtain ⟨s, hs⟩ := exists_isSStar hc hw
  exact ⟨s, hs, fun t ht => isSStar_unique hc hw ht hs⟩

theorem isSStar_sStar {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    IsSStar c w (sStar c w) := by
  unfold sStar
  rw [dif_pos (existsUnique_isSStar hc hw)]
  exact (existsUnique_isSStar hc hw).choose_spec.1

theorem sStar_mem {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    sStar c w ∈ Set.Ioo (sLo w) 0 := (isSStar_sStar hc hw).1

theorem phi_sStar {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    phi c w (sStar c w) = 0 := (isSStar_sStar hc hw).2

theorem sStar_eq {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {s : ℝ}
    (hs : IsSStar c w s) : sStar c w = s :=
  isSStar_unique hc hw (isSStar_sStar hc hw) hs

theorem phi_pos_of_lt {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {s : ℝ}
    (hs : sStar c w < s) (hs0 : s < 0) : 0 < phi c w s := by
  have h := phi_strictMonoOn hc hw (sStar_mem hc hw) ⟨lt_trans (sStar_mem hc hw).1 hs, hs0⟩ hs
  rwa [phi_sStar hc hw] at h

theorem phi_neg_of_lt {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {s : ℝ}
    (hs : sLo w < s) (hs0 : s < sStar c w) : phi c w s < 0 := by
  have h := phi_strictMonoOn hc hw ⟨hs, lt_trans hs0 (sStar_mem hc hw).2⟩ (sStar_mem hc hw) hs0
  rwa [phi_sStar hc hw] at h

/-! ### `zfun`: derivative, monotonicity on both sides of `sStar`, the edge `bHet` -/

theorem phi_eq_mul (c w : Fin M → ℝ) {s : ℝ} (hs0 : s ≠ 0) :
    phi c w s = s ^ 2 * zfunDeriv c w s := by
  unfold phi zfunDeriv
  rw [mul_sub, Finset.mul_sum]
  congr 1
  · field_simp
  · exact Finset.sum_congr rfl fun i _ => by ring

theorem zfunDeriv_pos {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {s : ℝ}
    (hs : sStar c w < s) (hs0 : s < 0) : 0 < zfunDeriv c w s := by
  have h := phi_pos_of_lt hc hw hs hs0
  rw [phi_eq_mul c w hs0.ne] at h
  exact (mul_pos_iff_of_pos_left (by nlinarith)).mp h

theorem zfunDeriv_neg {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {s : ℝ}
    (hs : sLo w < s) (hs0 : s < sStar c w) : zfunDeriv c w s < 0 := by
  have h := phi_neg_of_lt hc hw hs hs0
  have hs0' : s ≠ 0 := (lt_trans hs0 (sStar_mem hc hw).2).ne
  rw [phi_eq_mul c w hs0'] at h
  by_contra hcon
  have hcon' : 0 ≤ zfunDeriv c w s := not_lt.mp hcon
  nlinarith [mul_nonneg (sq_nonneg s) hcon']

theorem zfunDeriv_sStar {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    zfunDeriv c w (sStar c w) = 0 := by
  have h := phi_eq_mul c w (sStar_mem hc hw).2.ne
  rw [phi_sStar hc hw] at h
  have : sStar c w ^ 2 ≠ 0 := pow_ne_zero 2 (sStar_mem hc hw).2.ne
  exact (mul_eq_zero.mp h.symm).resolve_left this

theorem hasDerivAt_zfun {c w : Fin M → ℝ} (hw : ∃ i, w i ≠ 0) {s : ℝ} (hs : sLo w < s)
    (hs0 : s ≠ 0) : HasDerivAt (zfun c w) (zfunDeriv c w s) s := by
  have h1 : HasDerivAt (fun y : ℝ => -1 / y) (1 / s ^ 2) s := by
    have := (hasDerivAt_const s (-1 : ℝ)).div (hasDerivAt_id' s) hs0
    exact this.congr_deriv (by ring)
  have h2 : ∀ i ∈ (Finset.univ : Finset (Fin M)),
      HasDerivAt (fun y : ℝ => c i * w i ^ 2 / (1 + w i ^ 2 * y))
        (-(c i * w i ^ 4 / (1 + w i ^ 2 * s) ^ 2)) s := by
    intro i _
    have hne := (one_add_mul_pos hw hs i).ne'
    have := (hasDerivAt_const s (c i * w i ^ 2)).div
      (((hasDerivAt_id' s).const_mul (w i ^ 2)).const_add 1) hne
    exact this.congr_deriv (by ring)
  unfold zfun zfunDeriv
  refine (h1.add (HasDerivAt.fun_sum h2)).congr_deriv ?_
  rw [Finset.sum_neg_distrib]
  ring

theorem zfun_strictMonoOn {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    StrictMonoOn (zfun c w) (Set.Ico (sStar c w) 0) := by
  have hlo := (sStar_mem hc hw).1
  refine strictMonoOn_of_deriv_pos (convex_Ico _ _) ?_ ?_
  · intro s hs
    exact (hasDerivAt_zfun hw (lt_of_lt_of_le hlo hs.1) hs.2.ne).continuousAt.continuousWithinAt
  · intro s hs
    rw [interior_Ico] at hs
    rw [(hasDerivAt_zfun hw (lt_trans hlo hs.1) hs.2.ne).deriv]
    exact zfunDeriv_pos hc hw hs.1 hs.2

theorem zfun_strictAntiOn {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    StrictAntiOn (zfun c w) (Set.Ioc (sLo w) (sStar c w)) := by
  have hst := sStar_mem hc hw
  refine strictAntiOn_of_deriv_neg (convex_Ioc _ _) ?_ ?_
  · intro s hs
    exact (hasDerivAt_zfun hw hs.1 (lt_of_le_of_lt hs.2 hst.2).ne).continuousAt.continuousWithinAt
  · intro s hs
    rw [interior_Ioc] at hs
    rw [(hasDerivAt_zfun hw hs.1 (lt_trans hs.2 hst.2).ne).deriv]
    exact zfunDeriv_neg hc hw hs.1 hs.2

/-- `zfun s ≥ -1/s` on the domain: the sum is nonnegative. -/
theorem neg_inv_le_zfun {c w : Fin M → ℝ} (hc : ∀ i, 0 ≤ c i) (hw : ∃ i, w i ≠ 0) {s : ℝ}
    (hs : sLo w < s) : -1 / s ≤ zfun c w s := by
  unfold zfun
  have : 0 ≤ ∑ i, c i * w i ^ 2 / (1 + w i ^ 2 * s) := Finset.sum_nonneg fun i _ =>
    div_nonneg (mul_nonneg (hc i) (sq_nonneg _)) (one_add_mul_pos hw hs i).le
  linarith

theorem bHet_lt_zfun {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {s : ℝ}
    (hs : sStar c w < s) (hs0 : s < 0) : bHet c w < zfun c w s :=
  zfun_strictMonoOn hc hw ⟨le_rfl, (sStar_mem hc hw).2⟩ ⟨hs.le, hs0⟩ hs

/-- `bHet` is the minimum of `zfun` on `(sLo, 0)`. -/
theorem bHet_le_zfun {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {s : ℝ}
    (hs : s ∈ Set.Ioo (sLo w) 0) : bHet c w ≤ zfun c w s := by
  rcases lt_or_ge s (sStar c w) with h | h
  · exact (zfun_strictAntiOn hc hw ⟨hs.1, h.le⟩ ⟨(sStar_mem hc hw).1, le_rfl⟩ h).le
  · rcases eq_or_lt_of_le h with h' | h'
    · rw [← h']; exact le_rfl
    · exact (bHet_lt_zfun hc hw h' hs.2).le

theorem bHet_pos {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) : 0 < bHet c w := by
  have hst := sStar_mem hc hw
  have h1 : 0 < -1 / sStar c w := div_pos_of_neg_of_neg (by norm_num) hst.2
  exact lt_of_lt_of_le h1 (neg_inv_le_zfun (fun i => (hc i).le) hw hst.1)

/-! ### `sPhys`: existence, uniqueness, monotonicity -/

theorem exists_isSPhys {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ}
    (hz : bHet c w < z) : ∃ s, IsSPhys c w z s := by
  have hst := sStar_mem hc hw
  have hz0 : 0 < z := lt_trans (bHet_pos hc hw) hz
  set b := max (sStar c w) (-1 / z) / 2 with hb
  have hmax_neg : max (sStar c w) (-1 / z) < 0 :=
    max_lt hst.2 (div_neg_of_neg_of_pos (by norm_num) hz0)
  have hb_neg : b < 0 := by rw [hb]; linarith
  have hb_gt : sStar c w < b := by
    rw [hb]; have := le_max_left (sStar c w) (-1 / z); linarith
  have hb_gt' : -1 / z < b := by
    rw [hb]; have := le_max_right (sStar c w) (-1 / z); linarith
  have hzb : z < zfun c w b := by
    have h1 : z < -1 / b := by
      rw [lt_div_iff_of_neg hb_neg, mul_comm]
      exact (div_lt_iff₀ hz0).mp hb_gt'
    exact lt_of_lt_of_le h1 (neg_inv_le_zfun (fun i => (hc i).le) hw (lt_trans hst.1 hb_gt))
  have hcont : ContinuousOn (zfun c w) (Set.Icc (sStar c w) b) := fun s hs =>
    (hasDerivAt_zfun hw (lt_of_lt_of_le hst.1 hs.1)
      (lt_of_le_of_lt hs.2 hb_neg).ne).continuousAt.continuousWithinAt
  have hz' : zfun c w (sStar c w) < z := hz
  obtain ⟨s, hs, hfs⟩ := intermediate_value_Ioo hb_gt.le hcont ⟨hz', hzb⟩
  exact ⟨s, ⟨hs.1, lt_trans hs.2 hb_neg⟩, hfs⟩

theorem isSPhys_unique {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z s t : ℝ}
    (hs : IsSPhys c w z s) (ht : IsSPhys c w z t) : s = t :=
  (zfun_strictMonoOn hc hw).injOn ⟨hs.1.1.le, hs.1.2⟩ ⟨ht.1.1.le, ht.1.2⟩ (hs.2.trans ht.2.symm)

theorem existsUnique_isSPhys {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ}
    (hz : bHet c w < z) : ∃! s, IsSPhys c w z s := by
  obtain ⟨s, hs⟩ := exists_isSPhys hc hw hz
  exact ⟨s, hs, fun t ht => isSPhys_unique hc hw ht hs⟩

theorem isSPhys_sPhys {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ}
    (hz : bHet c w < z) : IsSPhys c w z (sPhys c w z) := by
  unfold sPhys
  rw [dif_pos (existsUnique_isSPhys hc hw hz)]
  exact (existsUnique_isSPhys hc hw hz).choose_spec.1

theorem sPhys_eq {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z s : ℝ}
    (hs : IsSPhys c w z s) : sPhys c w z = s := by
  have hz : bHet c w < z := by rw [← hs.2]; exact bHet_lt_zfun hc hw hs.1.1 hs.1.2
  exact isSPhys_unique hc hw (isSPhys_sPhys hc hw hz) hs

theorem sPhys_mem {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ}
    (hz : bHet c w < z) : sPhys c w z ∈ Set.Ioo (sStar c w) 0 := (isSPhys_sPhys hc hw hz).1

theorem zfun_sPhys {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ}
    (hz : bHet c w < z) : zfun c w (sPhys c w z) = z := (isSPhys_sPhys hc hw hz).2

theorem sStar_lt_sPhys {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ}
    (hz : bHet c w < z) : sStar c w < sPhys c w z := (sPhys_mem hc hw hz).1

theorem sPhys_neg {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ}
    (hz : bHet c w < z) : sPhys c w z < 0 := (sPhys_mem hc hw hz).2

theorem sLo_lt_sPhys {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ}
    (hz : bHet c w < z) : sLo w < sPhys c w z :=
  lt_trans (sStar_mem hc hw).1 (sStar_lt_sPhys hc hw hz)

theorem one_add_mul_sPhys_pos {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    {z : ℝ} (hz : bHet c w < z) (i : Fin M) : 0 < 1 + w i ^ 2 * sPhys c w z :=
  one_add_mul_pos hw (sLo_lt_sPhys hc hw hz) i

theorem sPhys_strictMonoOn {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    StrictMonoOn (sPhys c w) (Set.Ioi (bHet c w)) := by
  intro z₁ hz₁ z₂ hz₂ hlt
  by_contra h
  have h' : sPhys c w z₂ ≤ sPhys c w z₁ := not_lt.mp h
  have := ((zfun_strictMonoOn hc hw).le_iff_le ⟨(sStar_lt_sPhys hc hw hz₂).le, sPhys_neg hc hw hz₂⟩
    ⟨(sStar_lt_sPhys hc hw hz₁).le, sPhys_neg hc hw hz₁⟩).mpr h'
  rw [zfun_sPhys hc hw hz₂, zfun_sPhys hc hw hz₁] at this
  exact absurd hlt (not_lt.mpr this)

theorem sPhys_image {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    Set.Ioo (sStar c w) 0 ⊆ sPhys c w '' Set.Ioi (bHet c w) := by
  intro s hs
  exact ⟨zfun c w s, bHet_lt_zfun hc hw hs.1 hs.2, sPhys_eq hc hw ⟨hs, rfl⟩⟩

theorem continuousAt_sPhys {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ}
    (hz : bHet c w < z) : ContinuousAt (sPhys c w) z :=
  continuousAt_of_monotoneOn_of_image_mem_nhds (sPhys_strictMonoOn hc hw).monotoneOn
    (Ioi_mem_nhds hz) (Filter.mem_of_superset
      (Ioo_mem_nhds (sStar_lt_sPhys hc hw hz) (sPhys_neg hc hw hz)) (sPhys_image hc hw))

/-- The derivative of the physical branch, `s'(z) = 1/zfun'(s(z))`. -/
theorem hasDerivAt_sPhys {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ}
    (hz : bHet c w < z) : HasDerivAt (sPhys c w) (sPhysDeriv c w z) z := by
  unfold sPhysDeriv
  refine HasDerivAt.of_local_left_inverse (continuousAt_sPhys hc hw hz)
    (hasDerivAt_zfun hw (sLo_lt_sPhys hc hw hz) (sPhys_neg hc hw hz).ne)
    (zfunDeriv_pos hc hw (sStar_lt_sPhys hc hw hz) (sPhys_neg hc hw hz)).ne' ?_
  filter_upwards [Ioi_mem_nhds hz] with y hy
  exact zfun_sPhys hc hw hy

theorem sPhysDeriv_pos {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ}
    (hz : bHet c w < z) : 0 < sPhysDeriv c w z :=
  inv_pos.mpr (zfunDeriv_pos hc hw (sStar_lt_sPhys hc hw hz) (sPhys_neg hc hw hz))

theorem continuousOn_zfunDeriv {c w : Fin M → ℝ} (hw : ∃ i, w i ≠ 0) :
    ContinuousOn (zfunDeriv c w) (Set.Ioo (sLo w) 0) := by
  unfold zfunDeriv
  refine ContinuousOn.sub ?_ (continuousOn_finsetSum _ fun i _ => ?_)
  · exact ContinuousOn.div continuousOn_const (by fun_prop) fun s hs => pow_ne_zero 2 hs.2.ne
  · exact ContinuousOn.div continuousOn_const (by fun_prop) fun s hs =>
      pow_ne_zero 2 (one_add_mul_pos hw hs.1 i).ne'

/-- `s(z) → s⋆` from the right as `z ↓ bHet`. -/
theorem tendsto_sPhys_bHet {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    Tendsto (sPhys c w) (𝓝[>] bHet c w) (𝓝[>] sStar c w) := by
  have hst := sStar_mem hc hw
  rw [tendsto_nhdsWithin_iff]
  refine ⟨?_, ?_⟩
  · rw [tendsto_order]
    refine ⟨fun a ha => ?_, fun b hb => ?_⟩
    · filter_upwards [self_mem_nhdsWithin] with z hz
      exact lt_trans ha (sStar_lt_sPhys hc hw hz)
    · set b' := min b (sStar c w / 2) with hb'
      have hb'1 : sStar c w < b' := lt_min hb (by linarith [hst.2])
      have hb'2 : b' < 0 := lt_of_le_of_lt (min_le_right _ _) (by linarith [hst.2])
      have hzb : bHet c w < zfun c w b' := bHet_lt_zfun hc hw hb'1 hb'2
      filter_upwards [Ioo_mem_nhdsGT hzb] with z hz
      have hlt : sPhys c w z < b' := by
        by_contra h
        have h' : b' ≤ sPhys c w z := not_lt.mp h
        have := ((zfun_strictMonoOn hc hw).le_iff_le ⟨hb'1.le, hb'2⟩
          ⟨(sStar_lt_sPhys hc hw hz.1).le, sPhys_neg hc hw hz.1⟩).mpr h'
        rw [zfun_sPhys hc hw hz.1] at this
        exact absurd hz.2 (not_lt.mpr this)
      exact lt_of_lt_of_le hlt (min_le_left _ _)
  · filter_upwards [self_mem_nhdsWithin] with z hz
    exact sStar_lt_sPhys hc hw hz

/-- The square-root edge: `s'(z) → +∞` as `z ↓ bHet` (the R3⁻ input). -/
theorem sPhysDeriv_tendsto_atTop {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    Tendsto (sPhysDeriv c w) (𝓝[>] bHet c w) atTop := by
  have hst := sStar_mem hc hw
  have h1 : Tendsto (fun z => zfunDeriv c w (sPhys c w z)) (𝓝[>] bHet c w) (𝓝[>] 0) := by
    rw [tendsto_nhdsWithin_iff]
    refine ⟨?_, ?_⟩
    · have hcont : ContinuousAt (zfunDeriv c w) (sStar c w) :=
        (continuousOn_zfunDeriv hw).continuousAt (Ioo_mem_nhds hst.1 hst.2)
      have := hcont.tendsto.comp (tendsto_nhdsWithin_iff.mp (tendsto_sPhys_bHet hc hw)).1
      rwa [zfunDeriv_sStar hc hw] at this
    · filter_upwards [self_mem_nhdsWithin] with z hz
      exact zfunDeriv_pos hc hw (sStar_lt_sPhys hc hw hz) (sPhys_neg hc hw hz)
  exact (tendsto_inv_nhdsGT_zero.comp h1).congr fun z => rfl

/-! ### The outlier `rhoHet` and the branch lemma (`eq:assumption4`, `main_paper.tex:1368`) -/

theorem exists_w_ne_zero_of_root {θ w : Fin M → ℝ} (h : ∃ g, IsGammaTop θ w g) :
    ∃ i, w i ≠ 0 := by
  obtain ⟨j, hj⟩ := exists_signal_of_root (gammaTop_of_exists h)
  exact ⟨j, right_ne_zero_of_mul hj⟩

theorem gammaTop_pos {θ w : Fin M → ℝ} (h : ∃ g, IsGammaTop θ w g) : 0 < gammaTop θ w :=
  lt_of_le_of_lt (wSqMax_nonneg w) (gammaTop_of_exists h).1

/-- `-1/γ₁ > sLo`, since `γ₁ > max_i w_i²`. -/
theorem sLo_lt_neg_inv_gammaTop {θ w : Fin M → ℝ} (h : ∃ g, IsGammaTop θ w g) :
    sLo w < -1 / gammaTop θ w := by
  have hW := wSqMax_pos (exists_w_ne_zero_of_root h)
  have := one_div_lt_one_div_of_lt hW (gammaTop_of_exists h).1
  unfold sLo
  rw [neg_div, neg_div]
  linarith

theorem neg_inv_gammaTop_neg {θ w : Fin M → ℝ} (h : ∃ g, IsGammaTop θ w g) :
    -1 / gammaTop θ w < 0 :=
  div_neg_of_neg_of_pos (by norm_num) (gammaTop_pos h)

/-- `ρ = zfun (-1/γ₁)`: the closed form of plan section 2.4. -/
theorem rhoHet_eq_zfun {θ c w : Fin M → ℝ} (h : ∃ g, IsGammaTop θ w g) :
    rhoHet θ c w = zfun c w (-1 / gammaTop θ w) := by
  have hne : gammaTop θ w ≠ 0 := (gammaTop_pos h).ne'
  unfold rhoHet zfun
  rw [mul_add, mul_one, Finset.mul_sum]
  congr 1
  · field_simp
  · refine Finset.sum_congr rfl fun i _ => ?_
    have hi := (gammaTop_sub_pos h i).ne'
    have hden : 1 + w i ^ 2 * (-1 / gammaTop θ w)
        = (gammaTop θ w - w i ^ 2) / gammaTop θ w := by
      field_simp; ring
    rw [hden, div_div_eq_mul_div]
    ring

/-- `phi (-1/γ₁) = η₁`: the summands become `c_i w_i⁴/(γ₁ - w_i²)²`. -/
theorem phi_neg_inv_gammaTop {θ c w : Fin M → ℝ} (h : ∃ g, IsGammaTop θ w g) :
    phi c w (-1 / gammaTop θ w) = eta1 θ c w := by
  have hne : gammaTop θ w ≠ 0 := (gammaTop_pos h).ne'
  unfold phi eta1
  congr 1
  refine Finset.sum_congr rfl fun i _ => ?_
  have hi := (gammaTop_sub_pos h i).ne'
  have hden : 1 + w i ^ 2 * (-1 / gammaTop θ w)
      = (gammaTop θ w - w i ^ 2) / gammaTop θ w := by
    field_simp; ring
  rw [hden, div_pow, div_pow, div_div_eq_mul_div]
  field_simp

/-- `zfun' (-1/γ₁) = γ₁² η₁`. -/
theorem zfunDeriv_neg_inv_gammaTop {θ c w : Fin M → ℝ} (h : ∃ g, IsGammaTop θ w g) :
    zfunDeriv c w (-1 / gammaTop θ w) = gammaTop θ w ^ 2 * eta1 θ c w := by
  have hne : gammaTop θ w ≠ 0 := (gammaTop_pos h).ne'
  have hne' : -1 / gammaTop θ w ≠ 0 := (neg_inv_gammaTop_neg h).ne
  have h1 := phi_eq_mul c w hne'
  rw [phi_neg_inv_gammaTop h] at h1
  rw [h1]
  field_simp

/-- The branch lemma: `eq:assumption4` holds exactly when `-1/γ₁` lies on the physical
branch `(s⋆, 0)`. With the junk `gammaTop = 0` the right side is `0 ∈ (s⋆, 0)`, false. -/
theorem assumption4_iff_branch {θ c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) :
    Assumption4 θ c w ↔ -1 / gammaTop θ w ∈ Set.Ioo (sStar c w) 0 := by
  constructor
  · rintro ⟨hex, hlt⟩
    have hw := exists_w_ne_zero_of_root hex
    have hmem : -1 / gammaTop θ w ∈ Set.Ioo (sLo w) 0 :=
      ⟨sLo_lt_neg_inv_gammaTop hex, neg_inv_gammaTop_neg hex⟩
    refine ⟨?_, hmem.2⟩
    have hphi : 0 < phi c w (-1 / gammaTop θ w) := by
      rw [phi_neg_inv_gammaTop hex]; exact (eta1_pos_iff θ c w).mpr hlt
    by_contra hcon
    have hle : -1 / gammaTop θ w ≤ sStar c w := not_lt.mp hcon
    have := ((phi_strictMonoOn hc hw).le_iff_le hmem (sStar_mem hc hw)).mpr hle
    rw [phi_sStar hc hw] at this
    linarith
  · rintro ⟨h1, h2⟩
    have hex : ∃ g, IsGammaTop θ w g := by
      by_contra hno
      rw [gammaTop_eq_zero hno] at h2
      simp at h2
    refine ⟨hex, ?_⟩
    have hw := exists_w_ne_zero_of_root hex
    have hmem : -1 / gammaTop θ w ∈ Set.Ioo (sLo w) 0 :=
      ⟨sLo_lt_neg_inv_gammaTop hex, neg_inv_gammaTop_neg hex⟩
    have := phi_strictMonoOn hc hw (sStar_mem hc hw) hmem h1
    rw [phi_sStar hc hw, phi_neg_inv_gammaTop hex] at this
    exact (eta1_pos_iff θ c w).mp this

theorem bHet_lt_rhoHet {θ c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (h4 : Assumption4 θ c w) :
    bHet c w < rhoHet θ c w := by
  have hb := (assumption4_iff_branch hc).mp h4
  rw [rhoHet_eq_zfun h4.1]
  exact bHet_lt_zfun hc (exists_w_ne_zero_of_root h4.1) hb.1 hb.2

theorem sPhys_rhoHet {θ c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (h4 : Assumption4 θ c w) :
    sPhys c w (rhoHet θ c w) = -1 / gammaTop θ w :=
  sPhys_eq hc (exists_w_ne_zero_of_root h4.1)
    ⟨(assumption4_iff_branch hc).mp h4, (rhoHet_eq_zfun h4.1).symm⟩

theorem gamHet_rhoHet {θ c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (h4 : Assumption4 θ c w) :
    gamHet c w (rhoHet θ c w) = gammaTop θ w := by
  unfold gamHet
  rw [sPhys_rhoHet hc h4]
  have := (gammaTop_pos h4.1).ne'
  field_simp

theorem sPhysDeriv_rhoHet {θ c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (h4 : Assumption4 θ c w) :
    sPhysDeriv c w (rhoHet θ c w) = 1 / (gammaTop θ w ^ 2 * eta1 θ c w) := by
  unfold sPhysDeriv
  rw [sPhys_rhoHet hc h4, zfunDeriv_neg_inv_gammaTop h4.1, one_div]

/-- `γ(z) > max_i w_i²` on the physical branch. -/
theorem wSqMax_lt_gamHet {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ}
    (hz : bHet c w < z) : wSqMax w < gamHet c w z := by
  have hW := wSqMax_pos hw
  have hs0 := sPhys_neg hc hw hz
  have h1 : -1 < sPhys c w z * wSqMax w := by
    have := sLo_lt_sPhys hc hw hz
    unfold sLo at this
    exact (div_lt_iff₀ hW).mp this
  unfold gamHet
  rw [lt_div_iff_of_neg hs0]
  linarith

/-! ### The secular function `1 + F` (plan section 2.4) -/

/-- `1 + F(z) = (γ/z) · secular θ w γ` with `γ = -1/s(z)`, for every `z > bHet`. -/
theorem one_add_F_eq {θ c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) {z : ℝ}
    (hz : bHet c w < z) :
    1 + Fhet θ c w z = gamHet c w z / z * secular θ w (gamHet c w z) := by
  have hz0 : 0 < z := lt_trans (bHet_pos hc hw) hz
  have hs0 : sPhys c w z < 0 := sPhys_neg hc hw hz
  have hp : ∀ i, 0 < 1 + w i ^ 2 * sPhys c w z := one_add_mul_sPhys_pos hc hw hz
  have hzf : z = zfun c w (sPhys c w z) := (zfun_sPhys hc hw hz).symm
  set s := sPhys c w z with hs
  have hγ : gamHet c w z = -1 / s := rfl
  have hg : ∀ i, z * ghet c w i z = -1 / (1 + w i ^ 2 * s) := fun i => by
    have := (hp i).ne'
    unfold ghet
    rw [← hs]
    field_simp
  have hA : z * (1 + Fhet θ c w z) = -1 / s - ∑ i, θ i ^ 2 * w i ^ 2 / (1 + w i ^ 2 * s) := by
    have e1 : z * Phihet θ c w z = -∑ i, θ i ^ 2 * w i ^ 2 / (1 + w i ^ 2 * s) := by
      unfold Phihet
      rw [Finset.mul_sum, ← Finset.sum_neg_distrib]
      exact Finset.sum_congr rfl fun i _ => by rw [mul_left_comm, hg i]; ring
    have e2 : z * Psihet c w z = -∑ i, c i * w i ^ 2 / (1 + w i ^ 2 * s) := by
      unfold Psihet
      rw [Finset.mul_sum, ← Finset.sum_neg_distrib]
      exact Finset.sum_congr rfl fun i _ => by rw [mul_left_comm, hg i]; ring
    have e3 : z = -1 / s + ∑ i, c i * w i ^ 2 / (1 + w i ^ 2 * s) := hzf
    unfold Fhet
    rw [mul_add, mul_add, mul_one, e1, e2]
    nth_rewrite 1 [e3]
    ring
  have hB : -1 / s * secular θ w (-1 / s)
      = -1 / s - ∑ i, θ i ^ 2 * w i ^ 2 / (1 + w i ^ 2 * s) := by
    unfold secular
    rw [mul_add, mul_one, Finset.mul_sum, sub_eq_add_neg, ← Finset.sum_neg_distrib]
    congr 1
    refine Finset.sum_congr rfl fun i _ => ?_
    have h1 := (hp i).ne'
    have h2 := hs0.ne
    have hden : w i ^ 2 - -1 / s = (1 + w i ^ 2 * s) / s := by field_simp; ring
    rw [hden, div_div_eq_mul_div]
    field_simp
  have hR : z * (gamHet c w z / z * secular θ w (gamHet c w z))
      = -1 / s * secular θ w (-1 / s) := by
    rw [hγ]
    field_simp
  apply mul_left_cancel₀ hz0.ne'
  rw [hR, hA, hB]

theorem one_add_F_rhoHet {θ c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (h4 : Assumption4 θ c w) :
    1 + Fhet θ c w (rhoHet θ c w) = 0 := by
  rw [one_add_F_eq hc (exists_w_ne_zero_of_root h4.1) (bHet_lt_rhoHet hc h4), gamHet_rhoHet hc h4,
    (gammaTop_of_exists h4.1).2, mul_zero]

/-- The outlier is the unique root of `1 + F` above the edge. -/
theorem one_add_F_eq_zero_iff {θ c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (h4 : Assumption4 θ c w)
    {z : ℝ} (hz : bHet c w < z) : 1 + Fhet θ c w z = 0 ↔ z = rhoHet θ c w := by
  have hw := exists_w_ne_zero_of_root h4.1
  constructor
  · intro h0
    rw [one_add_F_eq hc hw hz] at h0
    have hz0 : 0 < z := lt_trans (bHet_pos hc hw) hz
    have hγ : 0 < gamHet c w z := lt_of_le_of_lt (wSqMax_nonneg w) (wSqMax_lt_gamHet hc hw hz)
    have hsec : secular θ w (gamHet c w z) = 0 :=
      (mul_eq_zero.mp h0).resolve_left (div_pos hγ hz0).ne'
    have hroot : IsGammaTop θ w (gamHet c w z) := ⟨wSqMax_lt_gamHet hc hw hz, hsec⟩
    have hγeq : gamHet c w z = gammaTop θ w := (gammaTop_eq hroot).symm
    have hs : sPhys c w z = sPhys c w (rhoHet θ c w) := by
      rw [sPhys_rhoHet hc h4, ← hγeq]
      unfold gamHet
      have := (sPhys_neg hc hw hz).ne
      field_simp
    exact (sPhys_strictMonoOn hc hw).injOn hz (bHet_lt_rhoHet hc h4) hs
  · rintro rfl
    exact one_add_F_rhoHet hc h4

/-! ### The overlap identity `1/(ρ F'(ρ)) = L(w)` (`main_paper.tex:1433`) -/

/-- `ρ F'(ρ) = LwDen/η₁`, the sum identity of plan section 2.4 in the `Lean` closed forms. -/
theorem rhoHet_mul_FhetDeriv {θ c w : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (h4 : Assumption4 θ c w) :
    rhoHet θ c w * FhetDeriv θ c w (rhoHet θ c w) = LwDen θ w / eta1 θ c w := by
  have hex := h4.1
  have hw := exists_w_ne_zero_of_root hex
  have hγ : 0 < gammaTop θ w := gammaTop_pos hex
  have he : 0 < eta1 θ c w := (eta1_pos_iff θ c w).mpr h4.2
  have hρ : 0 < rhoHet θ c w := lt_trans (bHet_pos hc hw) (bHet_lt_rhoHet hc h4)
  have hD : ∀ i, 0 < gammaTop θ w - w i ^ 2 := gammaTop_sub_pos hex
  have hs : sPhys c w (rhoHet θ c w) = -1 / gammaTop θ w := sPhys_rhoHet hc h4
  have hs' : sPhysDeriv c w (rhoHet θ c w) = 1 / (gammaTop θ w ^ 2 * eta1 θ c w) :=
    sPhysDeriv_rhoHet hc h4
  set γ := gammaTop θ w with hγdef
  set e := eta1 θ c w with hedef
  set ρ := rhoHet θ c w with hρdef
  have hT : ∀ i, ρ * (θ i ^ 2 * w i ^ 2 * ghetDeriv c w i ρ)
        + ρ * (c i * w i ^ 2 * ghetDeriv c w i ρ)
      = γ / ρ * (θ i ^ 2 * w i ^ 2 / (γ - w i ^ 2) + c i * w i ^ 2 / (γ - w i ^ 2))
        + 1 / e * (θ i ^ 2 * w i ^ 4 / (γ - w i ^ 2) ^ 2 + c i * w i ^ 4 / (γ - w i ^ 2) ^ 2) := by
    intro i
    unfold ghetDeriv
    rw [hs, hs']
    have h1 := (hD i).ne'
    have h2 := hγ.ne'
    have h3 := he.ne'
    have h4 := hρ.ne'
    have hden : 1 + w i ^ 2 * (-1 / γ) = (γ - w i ^ 2) / γ := by field_simp; ring
    rw [hden]
    field_simp
  have hsumθ : ∑ i, θ i ^ 2 * w i ^ 2 / (γ - w i ^ 2) = 1 := sum_pOf hex
  have hsumc : ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2 = 1 - e := by
    rw [hedef]; unfold eta1; rw [← hγdef]; ring
  have hρeq : ρ = γ * (1 + ∑ i, c i * w i ^ 2 / (γ - w i ^ 2)) := hρdef
  have hLw : LwDen θ w = 1 + ∑ i, θ i ^ 2 * w i ^ 4 / (γ - w i ^ 2) ^ 2 := by
    unfold LwDen
    rw [← hγdef, Finset.mul_sum, ← hsumθ, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    have h1 := (hD i).ne'
    field_simp
    ring
  have hB : 0 < 1 + ∑ i, c i * w i ^ 2 / (γ - w i ^ 2) := by
    rw [hρeq] at hρ
    exact (mul_pos_iff_of_pos_left hγ).mp hρ
  unfold FhetDeriv PhihetDeriv PsihetDeriv
  rw [mul_add, Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib,
    Finset.sum_congr rfl fun i _ => hT i, Finset.sum_add_distrib, ← Finset.mul_sum,
    ← Finset.mul_sum, Finset.sum_add_distrib, Finset.sum_add_distrib, hsumθ, hsumc, hLw]
  rw [hρeq]
  have h2 := hγ.ne'
  have h3 := he.ne'
  have h5 := hB.ne'
  field_simp
  ring

/-- The overlap identity: the R5-het limit `1/(ρ F'(ρ))` is the paper's `L(w)`. -/
theorem overlap_identity_het {θ c w : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (h4 : Assumption4 θ c w) :
    1 / (rhoHet θ c w * FhetDeriv θ c w (rhoHet θ c w)) = Lw θ c w := by
  rw [rhoHet_mul_FhetDeriv hc h4, Lw, if_pos h4, one_div_div]

/-! ### The Sudakov-Fernique edge dominates the exact edge (plan section 3.6) -/

/-- `bHet ≤ bSF`: evaluate `zfun` at `s₀ = -1/(ω(ω + σ))`, `ω = max_i |w_i|`,
`σ = √(∑ c_i w_i²)`, and bound each summand by the one with the largest weight. -/
theorem bHet_le_bSF {c w : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) :
    bHet c w ≤ bSF c w := by
  have hW := wSqMax_pos hw
  have hsum : 0 < ∑ i, c i * w i ^ 2 := by
    obtain ⟨k, hk⟩ := hw
    refine Finset.sum_pos' (fun i _ => ?_) ⟨k, Finset.mem_univ _, ?_⟩
    · have := (hc i).le; positivity
    · have := hc k; positivity
  set ω := Real.sqrt (wSqMax w) with hω
  set σ := Real.sqrt (∑ i, c i * w i ^ 2) with hσ
  have hω0 : 0 < ω := Real.sqrt_pos.mpr hW
  have hσ0 : 0 < σ := Real.sqrt_pos.mpr hsum
  have hω2 : ω ^ 2 = wSqMax w := Real.sq_sqrt hW.le
  have hσ2 : σ ^ 2 = ∑ i, c i * w i ^ 2 := Real.sq_sqrt hsum.le
  set s₀ := -1 / (ω * (ω + σ)) with hs₀def
  have hs₀neg : s₀ < 0 := div_neg_of_neg_of_pos (by norm_num) (by positivity)
  have hs₀lo : sLo w < s₀ := by
    unfold sLo
    rw [← hω2, hs₀def, neg_div, neg_div, neg_lt_neg_iff]
    exact one_div_lt_one_div_of_lt (by positivity) (by nlinarith)
  have hWs : 1 + wSqMax w * s₀ = σ / (ω + σ) := by
    rw [← hω2, hs₀def]; field_simp; ring
  have hWs0 : 0 < 1 + wSqMax w * s₀ := by rw [hWs]; positivity
  have h1 : bHet c w ≤ zfun c w s₀ := bHet_le_zfun hc hw ⟨hs₀lo, hs₀neg⟩
  have h2 : zfun c w s₀ ≤ -1 / s₀ + (∑ i, c i * w i ^ 2) / (1 + wSqMax w * s₀) := by
    unfold zfun
    rw [Finset.sum_div]
    refine add_le_add le_rfl (Finset.sum_le_sum fun i _ => ?_)
    have hle : 1 + wSqMax w * s₀ ≤ 1 + w i ^ 2 * s₀ := by
      have := mul_le_mul_of_nonpos_right (le_wSqMax w i) hs₀neg.le; linarith
    have hci := (hc i).le
    exact div_le_div_of_nonneg_left (by positivity) hWs0 hle
  have h3 : -1 / s₀ + (∑ i, c i * w i ^ 2) / (1 + wSqMax w * s₀) = (ω + σ) ^ 2 := by
    rw [hWs, ← hσ2, hs₀def]
    field_simp
  have h4 : bSF c w = (ω + σ) ^ 2 := rfl
  rw [h4]
  linarith

/-! ### The `M = 1`, `w = 1` reduction to `RMT/MP.lean` -/

section Single

variable {c₀ θ₀ : ℝ}

theorem wSqMax_one_single : wSqMax (1 : Fin 1 → ℝ) = 1 := by
  unfold wSqMax
  simp

theorem zfun_single (s : ℝ) : zfun (fun _ : Fin 1 => c₀) 1 s = -1 / s + c₀ / (1 + s) := by
  simp [zfun]

theorem phi_single (s : ℝ) : phi (fun _ : Fin 1 => c₀) 1 s = 1 - c₀ * s ^ 2 / (1 + s) ^ 2 := by
  simp [phi]

/-- At `M = 1`, `w = 1`: `s⋆ = -1/(1 + √c)`, the left end of the range of `MP.m`. -/
theorem sStar_single (hc : 0 < c₀) :
    sStar (fun _ : Fin 1 => c₀) 1 = -(1 / (1 + Real.sqrt c₀)) := by
  obtain ⟨r, hr0, rfl⟩ : ∃ r : ℝ, 0 < r ∧ c₀ = r ^ 2 :=
    ⟨Real.sqrt c₀, Real.sqrt_pos.mpr hc, (Real.sq_sqrt hc.le).symm⟩
  rw [Real.sqrt_sq hr0.le]
  have h1 : (1 + r) ≠ 0 := by positivity
  refine sStar_eq (fun _ => hc) ⟨0, one_ne_zero⟩ ⟨⟨?_, ?_⟩, ?_⟩
  · unfold sLo
    rw [wSqMax_one_single]
    have : 1 / (1 + r) < 1 := by rw [div_lt_one (by positivity)]; linarith
    linarith
  · have : 0 < 1 / (1 + r) := by positivity
    linarith
  · rw [phi_single]
    have hden : 1 + -(1 / (1 + r)) = r / (1 + r) := by field_simp; ring
    rw [hden]
    field_simp
    ring

/-- At `M = 1`, `w = 1`: `bHet = bulkEdge c = (1 + √c)²`. -/
theorem bHet_single (hc : 0 < c₀) : bHet (fun _ : Fin 1 => c₀) 1 = bulkEdge c₀ := by
  obtain ⟨r, hr0, rfl⟩ : ∃ r : ℝ, 0 < r ∧ c₀ = r ^ 2 :=
    ⟨Real.sqrt c₀, Real.sqrt_pos.mpr hc, (Real.sq_sqrt hc.le).symm⟩
  unfold bHet bulkEdge
  rw [sStar_single hc, zfun_single, Real.sqrt_sq hr0.le]
  have h1 : (1 + r) ≠ 0 := by positivity
  have hden : 1 + -(1 / (1 + r)) = r / (1 + r) := by field_simp; ring
  rw [hden]
  field_simp

/-- At `M = 1`, `w = 1`: the physical branch is the MP transform `MP.m`. -/
theorem sPhys_single (hc : 0 < c₀) {z : ℝ} (hz : bulkEdge c₀ < z) :
    sPhys (fun _ : Fin 1 => c₀) 1 z = MP.m c₀ z := by
  have hmem := MP.m_mem_Ioo hc hz
  have hm0 : MP.m c₀ z < 0 := hmem.2
  have hm1 : -1 < MP.m c₀ z := MP.neg_one_lt_m hc hz
  have hq := MP.m_quadratic hc hz
  refine sPhys_eq (fun _ => hc) ⟨0, one_ne_zero⟩ ⟨?_, ?_⟩
  · rw [sStar_single hc]; exact hmem
  · rw [zfun_single]
    have h1 : MP.m c₀ z ≠ 0 := hm0.ne
    have h2 : 1 + MP.m c₀ z ≠ 0 := by linarith
    rw [div_add_div _ _ h1 h2, div_eq_iff (mul_ne_zero h1 h2)]
    linear_combination (-1 : ℝ) * hq

/-- At `M = 1`, `w = 1`: the outlier is `rhoSq`. -/
theorem rhoHet_single (hc : 0 < c₀) (hθ : 0 < θ₀) (h : c₀ < θ₀ ^ 4) :
    rhoHet (fun _ : Fin 1 => θ₀) (fun _ => c₀) 1 = rhoSq θ₀ c₀ := by
  have hg : gammaTop (fun _ : Fin 1 => θ₀) 1 = 1 + θ₀ ^ 2 := by
    rw [gammaTop_one ⟨0, hθ.ne'⟩]
    simp
  unfold rhoHet
  rw [hg, MP.rhoSq_eq hc hθ h]
  simp only [Fin.sum_univ_one, Pi.one_apply, one_pow, mul_one]
  have hθne : θ₀ ≠ 0 := hθ.ne'
  field_simp
  ring

end Single

end MPhet
end StackedSVD
