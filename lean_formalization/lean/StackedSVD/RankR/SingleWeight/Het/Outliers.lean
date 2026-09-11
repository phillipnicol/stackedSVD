/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.SingleWeight.Het.Scalars
import StackedSVD.RankR.SingleWeight.Het.Forms
import StackedSVD.RankR.SingleWeight.Het.Frame
import StackedSVD.RankR.SingleWeight.Het.Count
import StackedSVD.RankR.Het.Outliers
import StackedSVD.RankR.Het.Edge
import StackedSVD.LinAlg.SpecWindow

/-!
# Track G, unit G3: the outliers of the single-weight Gram matrix

Unit G3 of `notes/archive/trackG_plan.md` (2026-09-05). Mirrors `RankR/Het/Outliers.lean:387`
(`tendsto_measure_count_Ioi_tau_het`) and `:592`
(`tendstoInProb_normSq_specProj_Ioi_tau_het`) at general `R_i`, under
`SingleWeight.EigSep` instead of `Scalars.Assumption4`.

## The route in five lines

1. `EigSep` gives, for every `l`, an outlier `ρ_l = swRho c w (γ l)` above the bulk edge, a
   positive `ν_l = swNu θ R c w (γ l) (z l)`, and `StrictAnti ρ`, hence a margin `g` that
   separates the `r` outliers from each other and from the edge.
2. The frame is the full matrix `Z = zMat z` of eigenvectors, not a column selection: every
   component is supercritical here, so the `exists_sum_equiv_split_pred` split of the aligned
   file disappears.
3. `ResolventLimitsSW` (unit G1) gives the `cform` and `cform2` limits of the columns of `Q`;
   bilinearity (`cform_mulVec_right`, `cform2_mulVec_mulVec`) turns them into the limits
   `swFmat θ R c w x *ᵥ u` and `u ⬝ᵥ (swFmat2 θ R c w x *ᵥ v)` of the frame combinations,
   whose values at `x = ρ_l`, `u = v = z l` are `-z l` and `ν l` (unit G0).
4. `Frame.count_of_split_of_frame_gap` (unit G2c) counts the eigenvalues above `τ`, and
   `OutliersR.align_detZ` (unit G2) reads the overlap of the scaled column `Q_k` off the
   same frame.
5. The `Fin s` sub-frame of the first `s` columns carries the overlap; `sum_castLE_eq_sum_filter`
   reindexes its sum onto the filter `(l : ℕ) < s` of the statement.

Every declaration is proved in full.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

/-! ### 1. Deterministic scalars from `EigSep` -/

namespace SingleWeight

variable {M r : ℕ} {rk : Fin M → ℕ}

section EigSepScalars

variable {θ : (i : Fin M) → Fin (rk i) → ℝ}
  {R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ} {w c : Fin M → ℝ}
  {γ : Fin r → ℝ} {z : Fin r → EuclideanSpace ℝ (Fin r)}

/-- Every root of `EigSep` sits above the bulk edge: `thresh` is exactly the hypothesis of
`SingleWeight.bHet_lt_swRho`. -/
theorem EigSep.bHet_lt_swRho (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hsep : EigSep θ R w c γ z) (l : Fin r) :
    MPhet.bHet c w < swRho c w (γ l) :=
  SingleWeight.bHet_lt_swRho hc hw (hsep.root l).1 (hsep.thresh l)

/-- Every `ν_l` of `EigSep` is positive. -/
theorem EigSep.swNu_pos (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hsep : EigSep θ R w c γ z) (l : Fin r) :
    0 < swNu θ R c w (γ l) (z l) :=
  SingleWeight.swNu_pos hc hw (hsep.root l).1 (hsep.thresh l) (hsep.eigvec l).1
    (hsep.eigvec l).2

/-- The outliers inherit the strict order of the roots: `swRho` is strictly increasing on the
admissible range and `γ` is strictly decreasing. -/
theorem EigSep.strictAnti_swRho (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hsep : EigSep θ R w c γ z) : StrictAnti fun l => swRho c w (γ l) := by
  intro k l hkl
  exact swRho_lt_swRho hc hw (hsep.root l).1 (hsep.thresh l) (hsep.root k).1 (hsep.thresh k)
    (hsep.sorted hkl)

/-- Distinct components carry distinct outliers. -/
theorem EigSep.injective_swRho (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hsep : EigSep θ R w c γ z) : Function.Injective fun l => swRho c w (γ l) :=
  (EigSep.strictAnti_swRho hc hw hsep).injective

/-- `F(ρ_l) z_l = -z_l`, the outlier equation of unit G0 read on `EigSep`. -/
theorem EigSep.swFmat_mulVec (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hsep : EigSep θ R w c γ z) (l : Fin r) :
    swFmat θ R c w (swRho c w (γ l)) *ᵥ WithLp.ofLp (z l) = -WithLp.ofLp (z l) :=
  swFmat_mulVec_eigvec hc hw (hsep.root l).1 (hsep.thresh l) (hsep.eigvec l).2

/-- `z_l ⬝ F'(ρ_l) z_l = ν_l`, the `ν` identity of unit G0 read on `EigSep`. -/
theorem EigSep.qform_swFmat2 (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hsep : EigSep θ R w c γ z) (l : Fin r) :
    WithLp.ofLp (z l) ⬝ᵥ (swFmat2 θ R c w (swRho c w (γ l)) *ᵥ WithLp.ofLp (z l))
      = swNu θ R c w (γ l) (z l) :=
  qform_swFmat2_eq_swNu hc hw (hsep.root l).1 (hsep.thresh l) (hsep.eigvec l).2

end EigSepScalars

/-- **The margin of a finite injective family above a floor.** One positive `g` bounds every
gap `ρ l - b` from below and separates every pair of distinct values. The `Fin r` family may
be empty; the constant `1` in the finite set keeps it nonempty. -/
theorem exists_margin_of_strictAnti {r : ℕ} {ρ : Fin r → ℝ} {b : ℝ}
    (hb : ∀ l, b < ρ l) (hinj : Function.Injective ρ) :
    ∃ g, 0 < g ∧ (∀ l, b + g ≤ ρ l) ∧ ∀ k l, k ≠ l → g ≤ |ρ k - ρ l| := by
  classical
  set T : Finset ℝ := insert (1 : ℝ)
    ((Finset.univ.image fun l : Fin r => ρ l - b) ∪
      ((Finset.univ.filter fun q : Fin r × Fin r => q.1 ≠ q.2).image
        fun q => |ρ q.1 - ρ q.2|)) with hTdef
  have hne : T.Nonempty := ⟨1, Finset.mem_insert_self _ _⟩
  have hmem1 : ∀ l : Fin r, ρ l - b ∈ T := by
    intro l
    rw [hTdef]
    exact Finset.mem_insert_of_mem
      (Finset.mem_union_left _ (Finset.mem_image.mpr ⟨l, Finset.mem_univ _, rfl⟩))
  have hmem2 : ∀ k l : Fin r, k ≠ l → |ρ k - ρ l| ∈ T := by
    intro k l hkl
    rw [hTdef]
    refine Finset.mem_insert_of_mem (Finset.mem_union_right _ (Finset.mem_image.mpr ?_))
    exact ⟨(k, l), Finset.mem_filter.mpr ⟨Finset.mem_univ _, hkl⟩, rfl⟩
  have hpos : ∀ x ∈ T, 0 < x := by
    intro x hx
    rw [hTdef, Finset.mem_insert] at hx
    rcases hx with rfl | hx
    · norm_num
    · rw [Finset.mem_union] at hx
      rcases hx with hx | hx
      · obtain ⟨l, -, rfl⟩ := Finset.mem_image.mp hx
        linarith [hb l]
      · obtain ⟨q, hq, rfl⟩ := Finset.mem_image.mp hx
        rw [Finset.mem_filter] at hq
        exact abs_pos.mpr (sub_ne_zero.mpr fun hcon => hq.2 (hinj hcon))
  refine ⟨T.min' hne, hpos _ (T.min'_mem hne), ?_, ?_⟩
  · intro l
    have h1 : T.min' hne ≤ ρ l - b := Finset.min'_le _ _ (hmem1 l)
    linarith
  · intro k l hkl
    exact Finset.min'_le _ _ (hmem2 k l hkl)

/-! ### 2. The eigenvector matrix and its column bound -/

/-- `Z`, the matrix whose column `k` is the eigenvector `z k`. -/
def zMat {r : ℕ} (z : Fin r → EuclideanSpace ℝ (Fin r)) : Matrix (Fin r) (Fin r) ℝ :=
  fun l k => WithLp.ofLp (z k) l

theorem zMat_apply {r : ℕ} (z : Fin r → EuclideanSpace ℝ (Fin r)) (l k : Fin r) :
    zMat z l k = WithLp.ofLp (z k) l := rfl

/-- Column `k` of `zMat z` is `z k`. -/
theorem zMat_col {r : ℕ} (z : Fin r → EuclideanSpace ℝ (Fin r)) (k : Fin r) :
    (fun l => zMat z l k) = WithLp.ofLp (z k) := rfl

/-- The `ℓ¹` bound `Zb` of the frame matrix: the sum of the absolute entries. -/
def zBound {r : ℕ} (z : Fin r → EuclideanSpace ℝ (Fin r)) : ℝ :=
  ∑ k, ∑ l, |zMat z l k|

theorem zBound_nonneg {r : ℕ} (z : Fin r → EuclideanSpace ℝ (Fin r)) : 0 ≤ zBound z :=
  Finset.sum_nonneg fun _ _ => Finset.sum_nonneg fun _ _ => abs_nonneg _

theorem sum_abs_zMat_le_zBound {r : ℕ} (z : Fin r → EuclideanSpace ℝ (Fin r)) (k : Fin r) :
    ∑ l, |zMat z l k| ≤ zBound z :=
  Finset.single_le_sum (f := fun k => ∑ l, |zMat z l k|)
    (fun _ _ => Finset.sum_nonneg fun _ _ => abs_nonneg _) (Finset.mem_univ k)

theorem abs_zMat_le_zBound {r : ℕ} (z : Fin r → EuclideanSpace ℝ (Fin r)) (l k : Fin r) :
    |zMat z l k| ≤ zBound z :=
  le_trans (Finset.single_le_sum (f := fun l => |zMat z l k|)
    (fun _ _ => abs_nonneg _) (Finset.mem_univ l)) (sum_abs_zMat_le_zBound z k)

/-! ### 3. The two scalar identities of the frame limits -/

/-- The limit of the first-order frame form, read as one entry of `swFmat *ᵥ u`. -/
theorem swFmat_mulVec_apply (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (c w : Fin M → ℝ) (x : ℝ)
    (u : Fin r → ℝ) (j : Fin r) :
    ∑ l, u l * ((∑ i, w i ^ 2 * MPhet.ghet c w i x * sigMat θ R i j l)
        + if j = l then MPhet.Psihet c w x else 0)
      = (swFmat θ R c w x *ᵥ u) j := by
  simp only [swFmat, Matrix.mulVec, dotProduct, Matrix.add_apply, Matrix.sum_apply,
    Matrix.smul_apply, Matrix.one_apply, smul_eq_mul, mul_ite, mul_one, mul_zero]
  exact Finset.sum_congr rfl fun l _ => by ring

/-- The limit of the second-order frame form, read as the quadratic form of `swFmat2`. -/
theorem swFmat2_qform_apply (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (c w : Fin M → ℝ) (x : ℝ)
    (u v : Fin r → ℝ) :
    ∑ j, ∑ l, u j * v l * ((∑ i, w i ^ 2 * MPhet.ghetDeriv c w i x * sigMat θ R i j l)
        + if j = l then MPhet.PsihetDeriv c w x else 0)
      = u ⬝ᵥ (swFmat2 θ R c w x *ᵥ v) := by
  simp only [swFmat2, Matrix.mulVec, dotProduct, Matrix.add_apply, Matrix.sum_apply,
    Matrix.smul_apply, Matrix.one_apply, smul_eq_mul, mul_ite, mul_one, mul_zero,
    Finset.mul_sum]
  exact Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun l _ => by ring

/-! ### 4. Reindexing the first `s` components -/

/-- The sum over the first `s` components, written on the filter of the statement. -/
theorem sum_castLE_eq_sum_filter {r s : ℕ} (hsr : s ≤ r) (f : Fin r → ℝ) :
    ∑ k : Fin s, f (Fin.castLE hsr k)
      = ∑ l ∈ Finset.univ.filter (fun l : Fin r => (l : ℕ) < s), f l := by
  classical
  have himg : (Finset.univ : Finset (Fin s)).image (Fin.castLE hsr)
      = Finset.univ.filter (fun l : Fin r => (l : ℕ) < s) := by
    ext l
    constructor
    · intro hl
      obtain ⟨k, -, rfl⟩ := Finset.mem_image.mp hl
      exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, k.isLt⟩
    · intro hl
      refine Finset.mem_image.mpr ⟨⟨(l : ℕ), (Finset.mem_filter.mp hl).2⟩,
        Finset.mem_univ _, ?_⟩
      rfl
  rw [← himg, Finset.sum_image fun a _ b _ hab => Fin.castLE_injective hsr hab]

/-- A nonnegative sum over the first `s` components is at most the full sum. -/
theorem sum_castLE_le {r s : ℕ} (hsr : s ≤ r) {f : Fin r → ℝ} (hf : ∀ l, 0 ≤ f l) :
    ∑ k : Fin s, f (Fin.castLE hsr k) ≤ ∑ l, f l := by
  classical
  rw [sum_castLE_eq_sum_filter hsr f]
  exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _) fun l _ _ => hf l

/-! ### 5. One accuracy for every bound -/

/-- **The accuracy `η` that drives the count and the overlap.** `CR` is the residual constant
`OutliersR.resCG`, `CG` the Gram constant `OutliersR.gramC`, `Zb` the column bound of the
frame, `mg` the margin above `τ`, `g` the margin between the outliers, and `ξ` any extra
accuracy the caller needs (`1` for the count, an `ε`-dependent bound for the overlap). The
five conclusions are exactly the smallness hypotheses of
`Frame.count_of_split_of_frame_gap` and of `OutliersR.align_detZ`. -/
theorem exists_accuracy_sw {r s : ℕ} {CR CG Zb mg g ξ : ℝ} (hsr : s ≤ r)
    (hCR : 0 ≤ CR) (hCG : 0 ≤ CG) (hZb0 : 0 ≤ Zb) (hmg : 0 < mg) (hg : 0 < g) (hξ : 0 < ξ) :
    ∃ η, 0 < η ∧ η ≤ ξ ∧ 2 * (CR * η) ≤ mg ∧ 4 * (CR * η) < g
      ∧ CG * (1 + Zb) * η ≤ 1 / 2
      ∧ (s : ℝ) * (CG * (1 + Zb) * η + (CR * η / mg) ^ 2) < 1
      ∧ (s : ℝ) * (CG * (1 + Zb) * η + CR * η / mg) ≤ 1 / 2 := by
  set A : ℝ := CG * (1 + Zb) with hAdef
  have hA0 : 0 ≤ A := mul_nonneg hCG (by linarith)
  set X : ℝ := (r : ℝ) * (A + CR / mg) with hXdef
  have hX0 : 0 ≤ X := mul_nonneg (Nat.cast_nonneg r) (add_nonneg hA0 (div_nonneg hCR hmg.le))
  have hCR1 : (0 : ℝ) < CR + 1 := by linarith
  set η : ℝ := min ξ (min (mg / (2 * (CR + 1))) (min (g / (8 * (CR + 1)))
      (min (1 / (2 * (A + 1))) (1 / (4 * (X + 1)))))) with hηdef
  have hη0 : 0 < η := by
    rw [hηdef]
    exact lt_min hξ (lt_min (div_pos hmg (by linarith))
      (lt_min (div_pos hg (by linarith))
        (lt_min (div_pos one_pos (by linarith)) (div_pos one_pos (by linarith)))))
  have hηξ : η ≤ ξ := by rw [hηdef]; exact min_le_left _ _
  have hη1 : η ≤ mg / (2 * (CR + 1)) := by
    rw [hηdef]; exact le_trans (min_le_right _ _) (min_le_left _ _)
  have hη2 : η ≤ g / (8 * (CR + 1)) := by
    rw [hηdef]
    exact le_trans (min_le_right _ _) (le_trans (min_le_right _ _) (min_le_left _ _))
  have hη3 : η ≤ 1 / (2 * (A + 1)) := by
    rw [hηdef]
    exact le_trans (min_le_right _ _) (le_trans (min_le_right _ _)
      (le_trans (min_le_right _ _) (min_le_left _ _)))
  have hη4 : η ≤ 1 / (4 * (X + 1)) := by
    rw [hηdef]
    exact le_trans (min_le_right _ _) (le_trans (min_le_right _ _)
      (le_trans (min_le_right _ _) (min_le_right _ _)))
  have hδmg : CR * η ≤ mg / 2 := by
    refine OutliersR.mul_le_of_le_div_add_one hCR ?_ (by linarith)
    rw [div_div]
    exact hη1
  have hδg : CR * η ≤ g / 8 := by
    refine OutliersR.mul_le_of_le_div_add_one hCR ?_ (by linarith)
    rw [div_div]
    exact hη2
  have hAη : A * η ≤ 1 / 2 :=
    le_trans (mul_le_mul_of_nonneg_left hη3 hA0) (OutliersR.mul_inv_le_half hA0)
  have hXη : X * η ≤ 1 / 4 := by
    have h1 : X * η ≤ X * (1 / (4 * (X + 1))) := mul_le_mul_of_nonneg_left hη4 hX0
    have h2 : X * (1 / (4 * (X + 1))) ≤ 1 / 4 := by
      rw [mul_one_div, div_le_div_iff₀ (by linarith) (by norm_num)]
      linarith
    linarith
  have hB0 : 0 ≤ CR * η / mg := div_nonneg (mul_nonneg hCR hη0.le) hmg.le
  have hB : CR * η / mg ≤ 1 / 2 := by
    rw [div_le_iff₀ hmg]
    linarith
  have hBsq : (CR * η / mg) ^ 2 ≤ CR * η / mg := by nlinarith
  have hs : (s : ℝ) ≤ (r : ℝ) := Nat.cast_le.mpr hsr
  have hs0 : (0 : ℝ) ≤ (s : ℝ) := Nat.cast_nonneg s
  have hAη0 : 0 ≤ A * η := mul_nonneg hA0 hη0.le
  have hkey : (r : ℝ) * (A * η + CR * η / mg) = X * η := by rw [hXdef]; ring
  have hmain : (s : ℝ) * (A * η + CR * η / mg) ≤ 1 / 4 := by
    have h1 : (s : ℝ) * (A * η + CR * η / mg) ≤ (r : ℝ) * (A * η + CR * η / mg) :=
      mul_le_mul_of_nonneg_right hs (by linarith)
    rw [hkey] at h1
    linarith
  refine ⟨η, hη0, hηξ, by linarith, by linarith, hAη, ?_, by linarith⟩
  have h1 : (s : ℝ) * (A * η + (CR * η / mg) ^ 2) ≤ (s : ℝ) * (A * η + CR * η / mg) :=
    mul_le_mul_of_nonneg_left (by linarith) hs0
  linarith

end SingleWeight

/-! ### 6. Bilinearity of the two resolvent forms -/

namespace OutliersR

open R4

variable {p r s : ℕ}

/-- The `k`-th column of `Q * Z` is `Q` applied to the `k`-th column of `Z`. Public copy of
the private `mulVec_col_eq` of `RankR/SingleWeight/Het/Frame.lean`. -/
theorem mul_col_eq_mulVec (Q : Matrix (Fin p) (Fin r) ℝ) (Z : Matrix (Fin r) (Fin s) ℝ)
    (k : Fin s) : (fun i => (Q * Z) i k) = Q *ᵥ fun l => Z l k := rfl

/-- Bilinearity of `cform` in its second argument along a column combination, the mirror of
the private `cform_mulVec_left` of `RankR/SingleWeight/Het/Frame.lean`. -/
theorem cform_mulVec_right (W : Matrix (Fin p) (Fin p) ℝ) (z : ℝ)
    (Q : Matrix (Fin p) (Fin r) ℝ) (x : Fin p → ℝ) (u : Fin r → ℝ) :
    cform W z x (Q *ᵥ u) = ∑ j, u j * cform W z x (fun i => Q i j) := by
  have h1 : cform W z x (Q *ᵥ u) = (x ᵥ* (resolv W z * Q)) ⬝ᵥ u := by
    change x ⬝ᵥ (resolv W z *ᵥ (Q *ᵥ u)) = (x ᵥ* (resolv W z * Q)) ⬝ᵥ u
    rw [Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec]
  have h2 : x ᵥ* (resolv W z * Q) = fun j => cform W z x (fun i => Q i j) := rfl
  rw [h1, h2]
  exact Finset.sum_congr rfl fun j _ => mul_comm _ _

/-- Bilinearity of `cform2` in its first argument along a column combination. -/
private theorem cform2_mulVec_left (W : Matrix (Fin p) (Fin p) ℝ) (z : ℝ)
    (Q : Matrix (Fin p) (Fin r) ℝ) (u : Fin r → ℝ) (y : Fin p → ℝ) :
    cform2 W z (Q *ᵥ u) y = ∑ j, u j * cform2 W z (fun i => Q i j) y := by
  have h1 : cform2 W z (Q *ᵥ u) y
      = (((resolv W z * resolv W z) *ᵥ y) ᵥ* Q) ⬝ᵥ u :=
    (dotProduct_comm (Q *ᵥ u) ((resolv W z * resolv W z) *ᵥ y)).trans
      (Matrix.dotProduct_mulVec _ _ _)
  have h2 : ((resolv W z * resolv W z) *ᵥ y) ᵥ* Q
      = fun j => cform2 W z (fun i => Q i j) y :=
    funext fun j => dotProduct_comm _ _
  rw [h1, h2]
  exact Finset.sum_congr rfl fun j _ => mul_comm _ _

/-- Bilinearity of `cform2` in its second argument along a column combination. -/
private theorem cform2_mulVec_right (W : Matrix (Fin p) (Fin p) ℝ) (z : ℝ)
    (Q : Matrix (Fin p) (Fin r) ℝ) (x : Fin p → ℝ) (v : Fin r → ℝ) :
    cform2 W z x (Q *ᵥ v) = ∑ l, v l * cform2 W z x (fun i => Q i l) := by
  have h1 : cform2 W z x (Q *ᵥ v) = (x ᵥ* (resolv W z * resolv W z * Q)) ⬝ᵥ v := by
    change x ⬝ᵥ ((resolv W z * resolv W z) *ᵥ (Q *ᵥ v))
      = (x ᵥ* (resolv W z * resolv W z * Q)) ⬝ᵥ v
    rw [Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec]
  have h2 : x ᵥ* (resolv W z * resolv W z * Q)
      = fun l => cform2 W z x (fun i => Q i l) := rfl
  rw [h1, h2]
  exact Finset.sum_congr rfl fun l _ => mul_comm _ _

/-- The double expansion of `cform2` along two column combinations. -/
theorem cform2_mulVec_mulVec (W : Matrix (Fin p) (Fin p) ℝ) (z : ℝ)
    (Q : Matrix (Fin p) (Fin r) ℝ) (u v : Fin r → ℝ) :
    cform2 W z (Q *ᵥ u) (Q *ᵥ v)
      = ∑ j, ∑ l, u j * v l * cform2 W z (fun i => Q i j) (fun i => Q i l) := by
  rw [cform2_mulVec_left]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [cform2_mulVec_right, Finset.mul_sum]
  exact Finset.sum_congr rfl fun l _ => by ring

/-! ### 7. The frame constants of a sub-family -/

/-- The residual constant of the first `s` components is at most the full one. -/
theorem resCG_castLE_le {r s : ℕ} (hsr : s ≤ r) (Cq : ℝ) (rr : ℕ) (ν : Fin r → ℝ) :
    resCG Cq rr (fun k => ν (Fin.castLE hsr k)) ≤ resCG Cq rr ν := by
  refine mul_le_mul_of_nonneg_left ?_ (Real.sqrt_nonneg _)
  exact SingleWeight.sum_castLE_le hsr fun l => inv_nonneg.mpr (Real.sqrt_nonneg _)

/-- The Gram constant of the first `s` components is at most the full one. -/
theorem gramC_castLE_le {r s : ℕ} (hsr : s ≤ r) (ρ ν : Fin r → ℝ) :
    gramC (fun k => ρ (Fin.castLE hsr k)) (fun k => ν (Fin.castLE hsr k)) ≤ gramC ρ ν := by
  have hT : ∀ k l : Fin r,
      0 ≤ (if ρ k = ρ l then (1 : ℝ) else 2 / |ρ k - ρ l|) * (Real.sqrt (ν k * ν l))⁻¹ := by
    intro k l
    refine mul_nonneg ?_ (by positivity)
    split
    · norm_num
    · positivity
  simp only [gramC]
  calc ∑ k : Fin s, ∑ l : Fin s,
        (if ρ (Fin.castLE hsr k) = ρ (Fin.castLE hsr l) then (1 : ℝ)
          else 2 / |ρ (Fin.castLE hsr k) - ρ (Fin.castLE hsr l)|)
          * (Real.sqrt (ν (Fin.castLE hsr k) * ν (Fin.castLE hsr l)))⁻¹
      ≤ ∑ k : Fin s, ∑ l : Fin r,
          (if ρ (Fin.castLE hsr k) = ρ l then (1 : ℝ) else 2 / |ρ (Fin.castLE hsr k) - ρ l|)
            * (Real.sqrt (ν (Fin.castLE hsr k) * ν l))⁻¹ :=
        Finset.sum_le_sum fun k _ =>
          SingleWeight.sum_castLE_le hsr fun l => hT (Fin.castLE hsr k) l
    _ ≤ ∑ k : Fin r, ∑ l : Fin r,
          (if ρ k = ρ l then (1 : ℝ) else 2 / |ρ k - ρ l|)
            * (Real.sqrt (ν k * ν l))⁻¹ :=
        SingleWeight.sum_castLE_le hsr fun k => Finset.sum_nonneg fun l _ => hT k l

/-- The sum of the inverse `ν` over the first `s` components is at most the full one. -/
theorem sum_inv_castLE_le {r s : ℕ} (hsr : s ≤ r) {ν : Fin r → ℝ} (hν : ∀ l, 0 < ν l) :
    ∑ k : Fin s, (ν (Fin.castLE hsr k))⁻¹ ≤ ∑ l, (ν l)⁻¹ :=
  SingleWeight.sum_castLE_le hsr fun l => inv_nonneg.mpr (hν l).le

end OutliersR

/-! ### 8. The column Gram limit -/

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- `N_l := ∑ w_i² (S_i)_{ll} + ∑ w_i² c_i`, the limit of `Q_l ⬝ᵥ Q_l` at general `R_i`. The
aligned twin is `UnalignedModelR.colGramLimit` (`RankR/Het/Outliers.lean:102`), where
`(S_i)_{ll}` is `θ_{il}²`. -/
noncomputable def colGramLimitSW (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ)
    (l : Fin r) : ℝ :=
  ∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i l l + ∑ i, w i ^ 2 * c i

/-- `0 < N_l`: the signal part is nonnegative because `(S_i)_{ll} = ∑_j (R_i)_{lj}² θ_{ij}²`,
and the noise floor `∑ w_i² c_i` is positive. -/
theorem colGramLimitSW_pos (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0) (l : Fin r) : 0 < m.colGramLimitSW w c l := by
  have h1 : 0 ≤ ∑ i, w i ^ 2 * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i l l := by
    refine Finset.sum_nonneg fun i _ => mul_nonneg (sq_nonneg _) ?_
    rw [sigMat_apply_gen]
    refine Finset.sum_nonneg fun j _ => ?_
    have h2 : m.R i l j * (m.tbl i).θ j ^ 2 * m.R i l j
        = m.R i l j ^ 2 * (m.tbl i).θ j ^ 2 := by ring
    rw [h2]
    positivity
  have h3 := colGramFloor_pos (c := c) (w := w) hc hw
  unfold colGramLimitSW
  linarith

/-- `Q_l ⬝ᵥ Q_l → N_l` in probability, the diagonal case of
`tendstoInProb_dotProduct_QmatHetR_col_gen` (unit G1). -/
theorem tendstoInProb_colGramSW [NeZero M] (m : UnalignedModelR μ M n d r rk)
    (w c : Fin M → ℝ) (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r) (l : Fin r) :
    TendstoInProb μ (fun N ω =>
        (fun q => m.QmatHetR w N ω q l) ⬝ᵥ (fun q => m.QmatHetR w N ω q l))
      (m.colGramLimitSW w c l) := by
  refine FormsR.tendstoInProb_congr_limit ?_
    (m.tendstoInProb_dotProduct_QmatHetR_col_gen w c hw hreg hG hpd l l)
  simp [colGramLimitSW]

/-! ### 9. The limits of the frame forms -/

/-- The limit of `Q_j ⬝ G₀(x) (Q u)` is the entry `j` of `swFmat θ R c w x *ᵥ u`: bilinearity
turns the column limits of `ResolventLimitsSW.cform_qcol` into the matrix product. -/
theorem tendstoInProb_cform_qcol_mulVec (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ)
    {b : ℝ}
    (H : m.ResolventLimitsSW w b
      (fun k l z => ∑ i, w i ^ 2 * MPhet.ghet c w i z
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
      (fun k l z => ∑ i, w i ^ 2 * MPhet.ghetDeriv c w i z
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
      (MPhet.Psihet c w) (MPhet.PsihetDeriv c w))
    (j : Fin r) (u : Fin r → ℝ) {x : ℝ} (hx : b < x) :
    TendstoInProb μ (fun N ω => R4.cform (m.W0hetR w N ω) x
        (fun q => m.QmatHetR w N ω q j) (m.QmatHetR w N ω *ᵥ u))
      ((SingleWeight.swFmat (fun i => (m.tbl i).θ) m.R c w x *ᵥ u) j) := by
  have hfun : (fun N (ω : Ω N) => R4.cform (m.W0hetR w N ω) x
        (fun q => m.QmatHetR w N ω q j) (m.QmatHetR w N ω *ᵥ u))
      = fun N ω => ∑ l, u l * R4.cform (m.W0hetR w N ω) x
          (fun q => m.QmatHetR w N ω q j) (fun q => m.QmatHetR w N ω q l) :=
    funext fun N => funext fun ω =>
      OutliersR.cform_mulVec_right (m.W0hetR w N ω) x (m.QmatHetR w N ω) _ u
  rw [hfun]
  have hlim := TendstoInProb.finsum (μ := μ)
    (F := fun l => fun N (ω : Ω N) => u l * R4.cform (m.W0hetR w N ω) x
      (fun q => m.QmatHetR w N ω q j) (fun q => m.QmatHetR w N ω q l))
    (v := fun l => u l * ((∑ i, w i ^ 2 * MPhet.ghet c w i x
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i j l)
      + if j = l then MPhet.Psihet c w x else 0))
    fun l => TendstoInProb.const_mul (u l) (H.cform_qcol j l hx)
  exact FormsR.tendstoInProb_congr_limit
    (SingleWeight.swFmat_mulVec_apply (fun i => (m.tbl i).θ) m.R c w x u j) hlim

/-- The limit of `(Q u) ⬝ G₀(x)² (Q v)` is the quadratic form of `swFmat2 θ R c w x`. -/
theorem tendstoInProb_cform2_qcol_mulVec (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ)
    {b : ℝ}
    (H : m.ResolventLimitsSW w b
      (fun k l z => ∑ i, w i ^ 2 * MPhet.ghet c w i z
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
      (fun k l z => ∑ i, w i ^ 2 * MPhet.ghetDeriv c w i z
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
      (MPhet.Psihet c w) (MPhet.PsihetDeriv c w))
    (u v : Fin r → ℝ) {x : ℝ} (hx : b < x) :
    TendstoInProb μ (fun N ω => R4.cform2 (m.W0hetR w N ω) x
        (m.QmatHetR w N ω *ᵥ u) (m.QmatHetR w N ω *ᵥ v))
      (u ⬝ᵥ (SingleWeight.swFmat2 (fun i => (m.tbl i).θ) m.R c w x *ᵥ v)) := by
  have hfun : (fun N (ω : Ω N) => R4.cform2 (m.W0hetR w N ω) x
        (m.QmatHetR w N ω *ᵥ u) (m.QmatHetR w N ω *ᵥ v))
      = fun N ω => ∑ j, ∑ l, u j * v l * R4.cform2 (m.W0hetR w N ω) x
          (fun q => m.QmatHetR w N ω q j) (fun q => m.QmatHetR w N ω q l) :=
    funext fun N => funext fun ω =>
      OutliersR.cform2_mulVec_mulVec (m.W0hetR w N ω) x (m.QmatHetR w N ω) u v
  rw [hfun]
  have hlim := TendstoInProb.finsum (μ := μ)
    (F := fun j => fun N (ω : Ω N) => ∑ l, u j * v l * R4.cform2 (m.W0hetR w N ω) x
      (fun q => m.QmatHetR w N ω q j) (fun q => m.QmatHetR w N ω q l))
    (v := fun j => ∑ l, u j * v l * ((∑ i, w i ^ 2 * MPhet.ghetDeriv c w i x
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i j l)
      + if j = l then MPhet.PsihetDeriv c w x else 0))
    fun j => TendstoInProb.finsum (μ := μ)
      (F := fun l => fun N (ω : Ω N) => u j * v l * R4.cform2 (m.W0hetR w N ω) x
        (fun q => m.QmatHetR w N ω q j) (fun q => m.QmatHetR w N ω q l))
      fun l => TendstoInProb.const_mul (u j * v l) (H.cform2_qcol j l hx)
  exact FormsR.tendstoInProb_congr_limit
    (SingleWeight.swFmat2_qform_apply (fun i => (m.tbl i).θ) m.R c w x u v) hlim

section Frame

variable {m : UnalignedModelR μ M n d r rk} {w c : Fin M → ℝ}
  {γ : Fin r → ℝ} {z : Fin r → EuclideanSpace ℝ (Fin r)}

/-- At `x = ρ_l` the first-order frame form tends to `-(z_l)_j`, by the outlier equation
`swFmat (ρ_l) *ᵥ z_l = -z_l` of unit G0. -/
theorem tendstoInProb_cform_qcol_frame
    (H : m.ResolventLimitsSW w (MPhet.bHet c w)
      (fun k l z => ∑ i, w i ^ 2 * MPhet.ghet c w i z
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
      (fun k l z => ∑ i, w i ^ 2 * MPhet.ghetDeriv c w i z
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
      (MPhet.Psihet c w) (MPhet.PsihetDeriv c w))
    (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z) (j l : Fin r) :
    TendstoInProb μ (fun N ω => R4.cform (m.W0hetR w N ω) (SingleWeight.swRho c w (γ l))
        (fun q => m.QmatHetR w N ω q j) (m.QmatHetR w N ω *ᵥ WithLp.ofLp (z l)))
      (-(WithLp.ofLp (z l) j)) := by
  refine FormsR.tendstoInProb_congr_limit ?_
    (m.tendstoInProb_cform_qcol_mulVec w c H j (WithLp.ofLp (z l))
      (SingleWeight.EigSep.bHet_lt_swRho hc hw hsep l))
  rw [SingleWeight.EigSep.swFmat_mulVec hc hw hsep l]
  rfl

/-- At `x = ρ_l` the second-order frame form tends to `ν_l`, by the `ν` identity of unit G0. -/
theorem tendstoInProb_cform2_qcol_frame
    (H : m.ResolventLimitsSW w (MPhet.bHet c w)
      (fun k l z => ∑ i, w i ^ 2 * MPhet.ghet c w i z
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
      (fun k l z => ∑ i, w i ^ 2 * MPhet.ghetDeriv c w i z
        * SingleWeight.sigMat (fun i => (m.tbl i).θ) m.R i k l)
      (MPhet.Psihet c w) (MPhet.PsihetDeriv c w))
    (hc : ∀ i, 0 < c i) (hw : ∃ i, w i ≠ 0)
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z) (l : Fin r) :
    TendstoInProb μ (fun N ω => R4.cform2 (m.W0hetR w N ω) (SingleWeight.swRho c w (γ l))
        (m.QmatHetR w N ω *ᵥ WithLp.ofLp (z l)) (m.QmatHetR w N ω *ᵥ WithLp.ofLp (z l)))
      (SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w (γ l) (z l)) :=
  FormsR.tendstoInProb_congr_limit (SingleWeight.EigSep.qform_swFmat2 hc hw hsep l)
    (m.tendstoInProb_cform2_qcol_mulVec w c H (WithLp.ofLp (z l)) (WithLp.ofLp (z l))
      (SingleWeight.EigSep.bHet_lt_swRho hc hw hsep l))

end Frame

/-! ### 10. The count of the eigenvalues above `τ` -/

/-- **The count of the eigenvalues above `τ`.** Under `EigSep` every component carries an
outlier `ρ_l = swRho c w (γ l)` above the bulk edge, and the outliers are strictly ordered, so
the sorted eigenvalues of `X_W X_Wᵀ` above a threshold `τ` that separates `ρ_{s-1}` from `ρ_s`
by the margin `mg` are exactly the indices below `s`, with probability tending to 1. The
aligned mirror is `tendsto_measure_count_Ioi_tau_het` (`RankR/Het/Outliers.lean:387`); the
`exists_sum_equiv_split_pred` split of that proof disappears, because `EigSep.thresh` makes
every component supercritical and the frame is the full matrix `zMat z`. -/
theorem tendsto_measure_count_Ioi_tau_sw [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)
    (hedge : m.HeteroEdgeR w c (MPhet.bHet c w))
    {γ : Fin r → ℝ} {z : Fin r → EuclideanSpace ℝ (Fin r)}
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z)
    {τ mg : ℝ} (hmg : 0 < mg) (hτ : MPhet.bHet c w < τ)
    {s : ℕ} (hsr : s ≤ r)
    (hgap : ∀ l : Fin r, ((l : ℕ) < s → τ + mg ≤ SingleWeight.swRho c w (γ l))
      ∧ (s ≤ (l : ℕ) → SingleWeight.swRho c w (γ l) + mg < τ)) :
    Tendsto (fun N => μ N {ω | ∀ k : Fin (Fintype.card (Fin (∑ i, n i N))),
      (τ < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ k
        ↔ (k : ℕ) < s)}) atTop (𝓝 1) := by
  have H := m.resolventLimitsSW_of_gaussian w c hc hw hreg hG hpd le_rfl hedge.edge
  -- 1. the frame scalars
  set ρv : Fin r → ℝ := fun l => SingleWeight.swRho c w (γ l) with hρvdef
  set νv : Fin r → ℝ :=
    fun l => SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w (γ l) (z l) with hνvdef
  set Z : Matrix (Fin r) (Fin r) ℝ := SingleWeight.zMat z with hZdef
  have hbρ : ∀ l, MPhet.bHet c w < ρv l := fun l =>
    SingleWeight.EigSep.bHet_lt_swRho hc hw hsep l
  have hνpos : ∀ l, 0 < νv l := fun l => SingleWeight.EigSep.swNu_pos hc hw hsep l
  have hρinj : Function.Injective ρv := SingleWeight.EigSep.injective_swRho hc hw hsep
  obtain ⟨g, hg0, hgb, hgsep⟩ := SingleWeight.exists_margin_of_strictAnti hbρ hρinj
  have habove : ∀ k : Fin r, (k : ℕ) < s → τ + mg ≤ ρv k := fun k hk => (hgap k).1 hk
  have hbelow0 : ∀ k : Fin r, s ≤ (k : ℕ) → ρv k + mg < τ := fun k hk => (hgap k).2 hk
  -- 2. the column-norm constant and the frame bound
  set Cq : ℝ := ∑ l, (m.colGramLimitSW w c l + 1) with hCqdef
  have hNl0 : ∀ l, 0 ≤ m.colGramLimitSW w c l + 1 := fun l =>
    (m.colGramLimitSW_pos w c hc hw l).le.trans (by linarith)
  have hCq0 : (0 : ℝ) ≤ Cq := Finset.sum_nonneg fun l _ => hNl0 l
  have hCqle : ∀ l, m.colGramLimitSW w c l + 1 ≤ Cq := fun l =>
    Finset.single_le_sum (f := fun l => m.colGramLimitSW w c l + 1) (fun l' _ => hNl0 l')
      (Finset.mem_univ l)
  set Zb : ℝ := SingleWeight.zBound z with hZbdef
  have hZb0 : (0 : ℝ) ≤ Zb := SingleWeight.zBound_nonneg z
  have hZbcol : ∀ k : Fin r, ∑ l, |Z l k| ≤ Zb := SingleWeight.sum_abs_zMat_le_zBound z
  -- 3. the accuracy
  have hCR0 : (0 : ℝ) ≤ OutliersR.resCG Cq r νv := OutliersR.resCG_nonneg Cq r νv
  have hCG0 : (0 : ℝ) ≤ OutliersR.gramC ρv νv := OutliersR.gramC_nonneg ρv νv
  obtain ⟨η, hη0, -, hδmg, hδg, hδ'half, hsmallc, -⟩ :=
    SingleWeight.exists_accuracy_sw (CR := OutliersR.resCG Cq r νv)
      (CG := OutliersR.gramC ρv νv) (Zb := Zb) (mg := mg) (g := g) (ξ := (1 : ℝ))
      hsr hCR0 hCG0 hZb0 hmg hg0 one_pos
  have hδ0 : (0 : ℝ) ≤ OutliersR.resCG Cq r νv * η := mul_nonneg hCR0 hη0.le
  have hδ'0 : (0 : ℝ) ≤ OutliersR.gramC ρv νv * (1 + Zb) * η :=
    mul_nonneg (mul_nonneg hCG0 (by linarith)) hη0.le
  -- 4. the four bad families
  set ε₀ : ℝ := min (g / 2) ((τ - MPhet.bHet c w) / 2) with hε₀def
  have hε₀0 : 0 < ε₀ := lt_min (by linarith) (by linarith)
  have hε₀g : ε₀ ≤ g / 2 := min_le_left _ _
  have hε₀τ : ε₀ ≤ (τ - MPhet.bHet c w) / 2 := min_le_right _ _
  have hedgeC : Tendsto (fun N => μ N {ω | lamMax (m.W0hetR w N ω)
      (m.isHermitian_W0hetR w N ω) ≤ MPhet.bHet c w + ε₀}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_lamMax_W0hetR_le w N
        (MPhet.bHet c w + ε₀)).nullMeasurableSet)
      (H.edge ε₀ hε₀0)
  have hcolT : ∀ l : Fin r, Tendsto (fun N => μ N
      {ω | (1 : ℝ) ≤ |(fun q => m.QmatHetR w N ω q l) ⬝ᵥ (fun q => m.QmatHetR w N ω q l)
        - m.colGramLimitSW w c l|}) atTop (𝓝 0) :=
    fun l => m.tendstoInProb_colGramSW w c hw hreg hG hpd l 1 one_pos
  have hE1T : ∀ q : Fin r × Fin r, Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform (m.W0hetR w N ω) (ρv q.2)
        (fun i => m.QmatHetR w N ω i q.1) (fun i => (m.QmatHetR w N ω * Z) i q.2)
        - -Z q.1 q.2|}) atTop (𝓝 0) := fun q =>
    tendstoInProb_cform_qcol_frame H hc hw hsep q.1 q.2 η hη0
  have hE2T : ∀ k : Fin r, Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform2 (m.W0hetR w N ω) (ρv k)
        (fun i => (m.QmatHetR w N ω * Z) i k) (fun i => (m.QmatHetR w N ω * Z) i k)
        - νv k|}) atTop (𝓝 0) := fun k =>
    tendstoInProb_cform2_qcol_frame H hc hw hsep k η hη0
  -- 5. assemble
  refine tendsto_measure_one_of_bad (s := fun N =>
      {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω)
          ≤ MPhet.bHet c w + ε₀}ᶜ
        ∪ ((⋃ l : Fin r, {ω | (1 : ℝ) ≤ |(fun q => m.QmatHetR w N ω q l)
                ⬝ᵥ (fun q => m.QmatHetR w N ω q l) - m.colGramLimitSW w c l|})
          ∪ ((⋃ q : Fin r × Fin r, {ω | η ≤ |R4.cform (m.W0hetR w N ω) (ρv q.2)
                (fun i => m.QmatHetR w N ω i q.1) (fun i => (m.QmatHetR w N ω * Z) i q.2)
                - -Z q.1 q.2|})
            ∪ (⋃ k : Fin r, {ω | η ≤ |R4.cform2 (m.W0hetR w N ω) (ρv k)
                (fun i => (m.QmatHetR w N ω * Z) i k) (fun i => (m.QmatHetR w N ω * Z) i k)
                - νv k|})))) ?_ ?_
  · intro N ω hω
    by_contra hbad
    simp only [Set.mem_union, Set.mem_compl_iff, Set.mem_iUnion, Set.mem_ofPred_eq,
      not_or, not_exists, not_le, not_lt] at hbad
    obtain ⟨hgood1, hgood2, hgood4, hgood5⟩ := hbad
    have hlamτ : lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) ≤ τ := by linarith
    have hcolb : ∀ l : Fin r, ∑ i, m.QmatHetR w N ω i l ^ 2 ≤ Cq := by
      intro l
      refine le_trans ?_ (hCqle l)
      have hsq : ∑ i, m.QmatHetR w N ω i l ^ 2
          = (fun q => m.QmatHetR w N ω q l) ⬝ᵥ (fun q => m.QmatHetR w N ω q l) :=
        Finset.sum_congr rfl fun i _ => sq _
      rw [hsq]
      have h1 := (abs_lt.mp (hgood2 l)).2
      linarith
    have hzedge : ∀ k : Fin r,
        lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω)
          + 2 * (OutliersR.resCG Cq r νv * η) < ρv k := by
      intro k
      have h1 := hgb k
      linarith
    have hzlt : ∀ k : Fin r,
        lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) < ρv k := by
      intro k
      have h1 := hzedge k
      linarith
    have hE1 : ∀ (l k : Fin r), |R4.cform (m.W0hetR w N ω) (ρv k)
        (fun i => m.QmatHetR w N ω i l) (fun i => (m.QmatHetR w N ω * Z) i k) + Z l k| ≤ η := by
      intro l k
      have h := (hgood4 (l, k)).le
      rwa [sub_neg_eq_add] at h
    have hE2 : ∀ k : Fin r, |R4.cform2 (m.W0hetR w N ω) (ρv k)
        (fun i => (m.QmatHetR w N ω * Z) i k) (fun i => (m.QmatHetR w N ω * Z) i k)
        - νv k| ≤ η := fun k => (hgood5 k).le
    have hres := fun k => OutliersR.residual_bound_Z (m.isHermitian_W0hetR w N ω)
      (m.gram_eq_hetR w N ω) hνpos hzlt hCq0 hη0.le hcolb hE1 k
    have hgram := OutliersR.gram_bound_Z (m.isHermitian_W0hetR w N ω) hρinj hνpos hzlt
      hZb0 hη0.le hZbcol hE1 hE2
    have hbelow : ∀ k : Fin r, s ≤ (k : ℕ) →
        ρv k + 2 * (OutliersR.resCG Cq r νv * η) ≤ τ := by
      intro k hk
      have h1 := hbelow0 k hk
      linarith
    have hsep' : ∀ k l : Fin r, s ≤ (k : ℕ) → s ≤ (l : ℕ) → k ≠ l →
        4 * (OutliersR.resCG Cq r νv * η) < |ρv k - ρv l| := by
      intro k l _ _ hkl
      have h1 := hgsep k l hkl
      linarith
    exact hω (Frame.count_of_split_of_frame_gap (m.isHermitian_W0hetR w N ω)
      (isHermitian_mul_transpose_self (m.stackXW w N ω)) (m.gram_eq_hetR w N ω) hsr
      hlamτ hmg hδ'0 hδ'half hzedge hbelow habove hsep' hres hgram hsmallc)
  · exact tendsto_measure_zero_union hedgeC
      (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hcolT)
        (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE1T)
          (tendsto_measure_zero_iUnion hE2T)))

/-! ### 11. The half-line overlap of a column of `Q` -/

/-- **The half-line overlap of a column of `Q`.** The squared norm of the projection of the
column `Q_k` on the eigenvalues of `X_W X_Wᵀ` above `τ` tends in probability to
`∑_{l < s} (z_l)_k² / ν_l`, the sum over the components whose outlier sits above `τ`. The
aligned mirror is `tendstoInProb_normSq_specProj_Ioi_tau_het`
(`RankR/Het/Outliers.lean:592`), where the sum has at most one term because the frame is a
selection of columns; here the frame is the matrix `zMat z` and the `s` surviving components
each contribute.

The route: `OutliersR.align_detZ` on the sub-frame of the first `s` columns of `zMat z`, with
the unit test vector `(√(Q_k ⬝ Q_k))⁻¹ • Q_k` on the event where the column Gram is close to
`colGramLimitSW w c k`; `OutliersR.normSq_specProj_eq_dot_mul` reads the answer back on `Q_k`,
and `SingleWeight.sum_castLE_eq_sum_filter` puts the sum on the filter of the statement. -/
theorem tendstoInProb_normSq_specProj_Ioi_tau_sw [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)
    (hedge : m.HeteroEdgeR w c (MPhet.bHet c w))
    {γ : Fin r → ℝ} {z : Fin r → EuclideanSpace ℝ (Fin r)}
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z)
    {τ mg : ℝ} (hmg : 0 < mg) (hτ : MPhet.bHet c w < τ)
    {s : ℕ} (hsr : s ≤ r)
    (hgap : ∀ l : Fin r, ((l : ℕ) < s → τ + mg ≤ SingleWeight.swRho c w (γ l))
      ∧ (s ≤ (l : ℕ) → SingleWeight.swRho c w (γ l) + mg < τ))
    (k : Fin r) :
    TendstoInProb μ
      (fun N ω => ‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ) (Set.Ioi τ)
        (WithLp.toLp 2 fun q => m.QmatHetR w N ω q k)‖ ^ 2)
      (∑ l ∈ Finset.univ.filter (fun l : Fin r => (l : ℕ) < s),
        (WithLp.ofLp (z l) k) ^ 2
          / SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w (γ l) (z l)) := by
  have H := m.resolventLimitsSW_of_gaussian w c hc hw hreg hG hpd le_rfl hedge.edge
  -- 1. the frame scalars, as in the count
  set ρv : Fin r → ℝ := fun l => SingleWeight.swRho c w (γ l) with hρvdef
  set νv : Fin r → ℝ :=
    fun l => SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w (γ l) (z l) with hνvdef
  set Z : Matrix (Fin r) (Fin r) ℝ := SingleWeight.zMat z with hZdef
  have hbρ : ∀ l, MPhet.bHet c w < ρv l := fun l =>
    SingleWeight.EigSep.bHet_lt_swRho hc hw hsep l
  have hνpos : ∀ l, 0 < νv l := fun l => SingleWeight.EigSep.swNu_pos hc hw hsep l
  have hρinj : Function.Injective ρv := SingleWeight.EigSep.injective_swRho hc hw hsep
  obtain ⟨g, hg0, hgb, hgsep⟩ := SingleWeight.exists_margin_of_strictAnti hbρ hρinj
  have habove : ∀ l : Fin r, (l : ℕ) < s → τ + mg ≤ ρv l := fun l hl => (hgap l).1 hl
  have hbelow0 : ∀ l : Fin r, s ≤ (l : ℕ) → ρv l + mg < τ := fun l hl => (hgap l).2 hl
  set Cq : ℝ := ∑ l, (m.colGramLimitSW w c l + 1) with hCqdef
  have hNl0 : ∀ l, 0 ≤ m.colGramLimitSW w c l + 1 := fun l =>
    (m.colGramLimitSW_pos w c hc hw l).le.trans (by linarith)
  have hCq0 : (0 : ℝ) ≤ Cq := Finset.sum_nonneg fun l _ => hNl0 l
  have hCqle : ∀ l, m.colGramLimitSW w c l + 1 ≤ Cq := fun l =>
    Finset.single_le_sum (f := fun l => m.colGramLimitSW w c l + 1) (fun l' _ => hNl0 l')
      (Finset.mem_univ l)
  set Zb : ℝ := SingleWeight.zBound z with hZbdef
  have hZb0 : (0 : ℝ) ≤ Zb := SingleWeight.zBound_nonneg z
  have hZbcol : ∀ l : Fin r, ∑ j, |Z j l| ≤ Zb := SingleWeight.sum_abs_zMat_le_zBound z
  have hCR0 : (0 : ℝ) ≤ OutliersR.resCG Cq r νv := OutliersR.resCG_nonneg Cq r νv
  have hCG0 : (0 : ℝ) ≤ OutliersR.gramC ρv νv := OutliersR.gramC_nonneg ρv νv
  -- 2. the sub-frame of the first `s` components
  set Zs : Matrix (Fin r) (Fin s) ℝ := Z.submatrix id (Fin.castLE hsr) with hZsdef
  set ρs : Fin s → ℝ := fun a => ρv (Fin.castLE hsr a) with hρsdef
  set νs : Fin s → ℝ := fun a => νv (Fin.castLE hsr a) with hνsdef
  have hρsinj : Function.Injective ρs := hρinj.comp (Fin.castLE_injective hsr)
  have hνspos : ∀ a, 0 < νs a := fun a => hνpos _
  have hρsτ : ∀ a : Fin s, τ + mg ≤ ρs a := fun a => habove _ a.isLt
  have hZbs : ∀ a : Fin s, ∑ j, |Zs j a| ≤ Zb := fun a => hZbcol _
  have hcmp1 : OutliersR.gramC ρs νs ≤ OutliersR.gramC ρv νv := by
    rw [hρsdef, hνsdef]
    exact OutliersR.gramC_castLE_le hsr ρv νv
  have hcmp2 : OutliersR.resCG Cq r νs ≤ OutliersR.resCG Cq r νv := by
    rw [hνsdef]
    exact OutliersR.resCG_castLE_le hsr Cq r νv
  have hcmp3 : ∑ a : Fin s, (νs a)⁻¹ ≤ ∑ l, (νv l)⁻¹ := by
    rw [hνsdef]
    exact OutliersR.sum_inv_castLE_le hsr hνpos
  -- 3. the scale of the column `k` and the frame scalars of the test vector
  have hNk : 0 < m.colGramLimitSW w c k := m.colGramLimitSW_pos w c hc hw k
  set sN : ℝ := Real.sqrt (m.colGramLimitSW w c k) with hsNdef
  have hsN : 0 < sN := Real.sqrt_pos.mpr hNk
  have hsN2 : sN ^ 2 = m.colGramLimitSW w c k := Real.sq_sqrt hNk.le
  have hsNi : 0 < sN⁻¹ := inv_pos.mpr hsN
  set Lb : ℝ := sN⁻¹ * Zb with hLbdef
  have hLb0 : (0 : ℝ) ≤ Lb := mul_nonneg hsNi.le hZb0
  set tv : Fin s → ℝ := fun a => sN⁻¹ * -Z k (Fin.castLE hsr a) with htvdef
  have hLbt : ∀ a, |tv a| ≤ Lb := by
    intro a
    have h1 : tv a = sN⁻¹ * -Z k (Fin.castLE hsr a) := rfl
    rw [h1, abs_mul, abs_of_pos hsNi, abs_neg, hLbdef]
    exact mul_le_mul_of_nonneg_left (SingleWeight.abs_zMat_le_zBound z k _) hsNi.le
  have htarget : m.colGramLimitSW w c k * ∑ a : Fin s, tv a ^ 2 / νs a
      = ∑ l ∈ Finset.univ.filter (fun l : Fin r => (l : ℕ) < s),
        (WithLp.ofLp (z l) k) ^ 2
          / SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w (γ l) (z l) := by
    have hterm : ∀ a : Fin s, tv a ^ 2 / νs a
        = (m.colGramLimitSW w c k)⁻¹
          * ((fun l : Fin r => (WithLp.ofLp (z l) k) ^ 2
              / SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w (γ l) (z l))
            (Fin.castLE hsr a)) := by
      intro a
      have h1 : tv a ^ 2 = (m.colGramLimitSW w c k)⁻¹ * Z k (Fin.castLE hsr a) ^ 2 := by
        have h2 : tv a = sN⁻¹ * -Z k (Fin.castLE hsr a) := rfl
        rw [h2, mul_pow, neg_sq, inv_pow, hsN2]
      rw [h1]
      exact mul_div_assoc _ _ _
    rw [Finset.sum_congr rfl fun a (_ : a ∈ Finset.univ) => hterm a, ← Finset.mul_sum,
      SingleWeight.sum_castLE_eq_sum_filter hsr
        (fun l : Fin r => (WithLp.ofLp (z l) k) ^ 2
          / SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w (γ l) (z l)),
      ← mul_assoc, mul_inv_cancel₀ hNk.ne', one_mul]
  -- 4. the accuracy, now `ε` dependent
  intro ε hε
  have hε2 : (0 : ℝ) < ε / 2 := by linarith
  have hνsum0 : (0 : ℝ) ≤ ∑ l, (νv l)⁻¹ :=
    Finset.sum_nonneg fun l _ => inv_nonneg.mpr (hνpos l).le
  set B : ℝ := 5 * (s : ℝ) * (OutliersR.gramC ρv νv * (1 + Zb)
      + OutliersR.resCG Cq r νv / mg) + (1 + 2 * Lb) * ∑ l, (νv l)⁻¹ with hBdef
  have hB0 : (0 : ℝ) ≤ B := by
    have h1 : (0 : ℝ) ≤ OutliersR.gramC ρv νv * (1 + Zb) + OutliersR.resCG Cq r νv / mg :=
      add_nonneg (mul_nonneg hCG0 (by linarith)) (div_nonneg hCR0 hmg.le)
    have h2 : (0 : ℝ) ≤ 5 * (s : ℝ) := by positivity
    have h3 : (0 : ℝ) ≤ (1 + 2 * Lb) * ∑ l, (νv l)⁻¹ :=
      mul_nonneg (by linarith) hνsum0
    rw [hBdef]
    nlinarith
  set K : ℝ := (m.colGramLimitSW w c k + 1) * B with hKdef
  have hK0 : (0 : ℝ) ≤ K := mul_nonneg (hNl0 k) hB0
  set Tv : ℝ := ∑ a : Fin s, tv a ^ 2 / νs a with hTvdef
  set η₂ : ℝ := min (m.colGramLimitSW w c k / 2) (ε / 2 / (2 * (|Tv| + 1))) with hη₂def
  have hη₂0 : 0 < η₂ := by
    rw [hη₂def]
    refine lt_min (by linarith) (div_pos hε2 ?_)
    have h1 := abs_nonneg Tv
    linarith
  have hη₂half : η₂ ≤ m.colGramLimitSW w c k / 2 := by rw [hη₂def]; exact min_le_left _ _
  have hη₂T : η₂ * |Tv| ≤ ε / 2 / 2 := by
    have h1 : η₂ ≤ ε / 2 / (2 * (|Tv| + 1)) := by rw [hη₂def]; exact min_le_right _ _
    rw [mul_comm]
    exact le_trans (mul_le_mul_of_nonneg_left h1 (abs_nonneg _))
      (OutliersR.mul_div_le_half (abs_nonneg _) hε2.le)
  obtain ⟨η, hη0, hηξ, hδmg, hδg, hδ'half, hsmallc, hsmalla⟩ :=
    SingleWeight.exists_accuracy_sw (CR := OutliersR.resCG Cq r νv)
      (CG := OutliersR.gramC ρv νv) (Zb := Zb) (mg := mg) (g := g)
      (ξ := min 1 (ε / 2 / (2 * (K + 1)))) hsr hCR0 hCG0 hZb0 hmg hg0
      (lt_min one_pos (div_pos hε2 (by linarith)))
  have hη1 : η ≤ 1 := le_trans hηξ (min_le_left _ _)
  have hηK : K * η ≤ ε / 2 / 2 :=
    le_trans (mul_le_mul_of_nonneg_left (le_trans hηξ (min_le_right _ _)) hK0)
      (OutliersR.mul_div_le_half hK0 hε2.le)
  have hδ0 : (0 : ℝ) ≤ OutliersR.resCG Cq r νv * η := mul_nonneg hCR0 hη0.le
  have hδ'0 : (0 : ℝ) ≤ OutliersR.gramC ρv νv * (1 + Zb) * η :=
    mul_nonneg (mul_nonneg hCG0 (by linarith)) hη0.le
  have hsmalls : (s : ℝ) * (OutliersR.gramC ρs νs * (1 + Zb) * η
      + OutliersR.resCG Cq r νs * η / mg) ≤ 1 / 2 := by
    refine le_trans ?_ hsmalla
    have h1 : OutliersR.gramC ρs νs * (1 + Zb) * η
        ≤ OutliersR.gramC ρv νv * (1 + Zb) * η :=
      mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right hcmp1 (by linarith)) hη0.le
    have h2 : OutliersR.resCG Cq r νs * η / mg ≤ OutliersR.resCG Cq r νv * η / mg := by
      rw [div_eq_mul_inv, div_eq_mul_inv]
      exact mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right hcmp2 hη0.le)
        (inv_pos.mpr hmg).le
    exact mul_le_mul_of_nonneg_left (by linarith) (Nat.cast_nonneg s)
  have hRHS : 5 * (s : ℝ) * (OutliersR.gramC ρs νs * (1 + Zb) * η
        + OutliersR.resCG Cq r νs * η / mg)
      + (1 + 2 * Lb) * (∑ a : Fin s, (νs a)⁻¹) * η ≤ B * η := by
    have hAs : OutliersR.gramC ρs νs * (1 + Zb) + OutliersR.resCG Cq r νs / mg
        ≤ OutliersR.gramC ρv νv * (1 + Zb) + OutliersR.resCG Cq r νv / mg := by
      have h1 : OutliersR.gramC ρs νs * (1 + Zb) ≤ OutliersR.gramC ρv νv * (1 + Zb) :=
        mul_le_mul_of_nonneg_right hcmp1 (by linarith)
      have h2 : OutliersR.resCG Cq r νs / mg ≤ OutliersR.resCG Cq r νv / mg := by
        rw [div_eq_mul_inv, div_eq_mul_inv]
        exact mul_le_mul_of_nonneg_right hcmp2 (inv_pos.mpr hmg).le
      linarith
    have hbase : 5 * (s : ℝ) * (OutliersR.gramC ρs νs * (1 + Zb)
          + OutliersR.resCG Cq r νs / mg) + (1 + 2 * Lb) * ∑ a : Fin s, (νs a)⁻¹ ≤ B := by
      have h1 : 5 * (s : ℝ) * (OutliersR.gramC ρs νs * (1 + Zb)
            + OutliersR.resCG Cq r νs / mg)
          ≤ 5 * (s : ℝ) * (OutliersR.gramC ρv νv * (1 + Zb)
            + OutliersR.resCG Cq r νv / mg) :=
        mul_le_mul_of_nonneg_left hAs (by positivity)
      have h2 : (1 + 2 * Lb) * (∑ a : Fin s, (νs a)⁻¹) ≤ (1 + 2 * Lb) * ∑ l, (νv l)⁻¹ :=
        mul_le_mul_of_nonneg_left hcmp3 (by linarith)
      rw [hBdef]
      linarith
    have he : 5 * (s : ℝ) * (OutliersR.gramC ρs νs * (1 + Zb) * η
          + OutliersR.resCG Cq r νs * η / mg)
        + (1 + 2 * Lb) * (∑ a : Fin s, (νs a)⁻¹) * η
        = (5 * (s : ℝ) * (OutliersR.gramC ρs νs * (1 + Zb)
            + OutliersR.resCG Cq r νs / mg)
          + (1 + 2 * Lb) * ∑ a : Fin s, (νs a)⁻¹) * η := by ring
    rw [he]
    exact mul_le_mul_of_nonneg_right hbase hη0.le
  -- 5. the six bad families
  set ε₀ : ℝ := min (g / 2) ((τ - MPhet.bHet c w) / 2) with hε₀def
  have hε₀0 : 0 < ε₀ := lt_min (by linarith) (by linarith)
  have hε₀g : ε₀ ≤ g / 2 := min_le_left _ _
  have hε₀τ : ε₀ ≤ (τ - MPhet.bHet c w) / 2 := min_le_right _ _
  have hedgeC : Tendsto (fun N => μ N {ω | lamMax (m.W0hetR w N ω)
      (m.isHermitian_W0hetR w N ω) ≤ MPhet.bHet c w + ε₀}ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_lamMax_W0hetR_le w N
        (MPhet.bHet c w + ε₀)).nullMeasurableSet)
      (H.edge ε₀ hε₀0)
  have hcolT : ∀ l : Fin r, Tendsto (fun N => μ N
      {ω | (1 : ℝ) ≤ |(fun q => m.QmatHetR w N ω q l) ⬝ᵥ (fun q => m.QmatHetR w N ω q l)
        - m.colGramLimitSW w c l|}) atTop (𝓝 0) :=
    fun l => m.tendstoInProb_colGramSW w c hw hreg hG hpd l 1 one_pos
  have hcol2T : Tendsto (fun N => μ N
      {ω | η₂ ≤ |(fun q => m.QmatHetR w N ω q k) ⬝ᵥ (fun q => m.QmatHetR w N ω q k)
        - m.colGramLimitSW w c k|}) atTop (𝓝 0) :=
    m.tendstoInProb_colGramSW w c hw hreg hG hpd k η₂ hη₂0
  have hE1T : ∀ q : Fin r × Fin r, Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform (m.W0hetR w N ω) (ρv q.2)
        (fun i => m.QmatHetR w N ω i q.1) (fun i => (m.QmatHetR w N ω * Z) i q.2)
        - -Z q.1 q.2|}) atTop (𝓝 0) := fun q =>
    tendstoInProb_cform_qcol_frame H hc hw hsep q.1 q.2 η hη0
  have hE2T : ∀ l : Fin r, Tendsto (fun N => μ N
      {ω | η ≤ |R4.cform2 (m.W0hetR w N ω) (ρv l)
        (fun i => (m.QmatHetR w N ω * Z) i l) (fun i => (m.QmatHetR w N ω * Z) i l)
        - νv l|}) atTop (𝓝 0) := fun l =>
    tendstoInProb_cform2_qcol_frame H hc hw hsep l η hη0
  have hE3T : ∀ a : Fin s, Tendsto (fun N => μ N
      {ω | η ≤ |(Real.sqrt ((fun q => m.QmatHetR w N ω q k)
            ⬝ᵥ (fun q => m.QmatHetR w N ω q k)))⁻¹
          * R4.cform (m.W0hetR w N ω) (ρs a) (fun q => m.QmatHetR w N ω q k)
            (fun i => (m.QmatHetR w N ω * Zs) i a)
        - tv a|}) atTop (𝓝 0) := by
    intro a
    have h2 := tendstoInProb_cform_qcol_frame H hc hw hsep k (Fin.castLE hsr a)
    have h3 : TendstoInProb μ (fun N (ω : Ω N) =>
        (Real.sqrt ((fun q => m.QmatHetR w N ω q k)
          ⬝ᵥ (fun q => m.QmatHetR w N ω q k)))⁻¹) sN⁻¹ :=
      (m.tendstoInProb_colGramSW w c hw hreg hG hpd k).comp_continuous
        (φ := fun x => (Real.sqrt x)⁻¹) (Real.continuous_sqrt.continuousAt.inv₀ hsN.ne')
    exact (h3.mul h2) η hη0
  -- 6. assemble
  refine tendsto_measure_zero_of_subset (t := fun N =>
      {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω)
          ≤ MPhet.bHet c w + ε₀}ᶜ
        ∪ ((⋃ l : Fin r, {ω | (1 : ℝ) ≤ |(fun q => m.QmatHetR w N ω q l)
                ⬝ᵥ (fun q => m.QmatHetR w N ω q l) - m.colGramLimitSW w c l|})
          ∪ ({ω | η₂ ≤ |(fun q => m.QmatHetR w N ω q k)
                ⬝ᵥ (fun q => m.QmatHetR w N ω q k) - m.colGramLimitSW w c k|}
            ∪ ((⋃ q : Fin r × Fin r, {ω | η ≤ |R4.cform (m.W0hetR w N ω) (ρv q.2)
                  (fun i => m.QmatHetR w N ω i q.1) (fun i => (m.QmatHetR w N ω * Z) i q.2)
                  - -Z q.1 q.2|})
              ∪ ((⋃ l : Fin r, {ω | η ≤ |R4.cform2 (m.W0hetR w N ω) (ρv l)
                    (fun i => (m.QmatHetR w N ω * Z) i l)
                    (fun i => (m.QmatHetR w N ω * Z) i l) - νv l|})
                ∪ (⋃ a : Fin s, {ω | η ≤ |(Real.sqrt ((fun q => m.QmatHetR w N ω q k)
                        ⬝ᵥ (fun q => m.QmatHetR w N ω q k)))⁻¹
                      * R4.cform (m.W0hetR w N ω) (ρs a) (fun q => m.QmatHetR w N ω q k)
                        (fun i => (m.QmatHetR w N ω * Zs) i a) - tv a|})))))) ?_ ?_
  · intro N ω hω
    by_contra hbad
    simp only [Set.mem_union, Set.mem_compl_iff, Set.mem_iUnion, Set.mem_ofPred_eq,
      not_or, not_exists, not_le, not_lt] at hbad
    obtain ⟨hgood1, hgood2, hgood3, hgood4, hgood5, hgood6⟩ := hbad
    have hlamτ : lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) ≤ τ := by linarith
    have hcolb : ∀ l : Fin r, ∑ i, m.QmatHetR w N ω i l ^ 2 ≤ Cq := by
      intro l
      refine le_trans ?_ (hCqle l)
      have hsq : ∑ i, m.QmatHetR w N ω i l ^ 2
          = (fun q => m.QmatHetR w N ω q l) ⬝ᵥ (fun q => m.QmatHetR w N ω q l) :=
        Finset.sum_congr rfl fun i _ => sq _
      rw [hsq]
      have h1 := (abs_lt.mp (hgood2 l)).2
      linarith
    have hgpos : 0 < (fun q => m.QmatHetR w N ω q k) ⬝ᵥ (fun q => m.QmatHetR w N ω q k) := by
      have h1 := (abs_lt.mp hgood3).1
      linarith
    have hgle : (fun q => m.QmatHetR w N ω q k) ⬝ᵥ (fun q => m.QmatHetR w N ω q k)
        ≤ m.colGramLimitSW w c k + 1 := by
      have h1 := (abs_lt.mp (hgood2 k)).2
      linarith
    have hzedge : ∀ l : Fin r,
        lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω)
          + 2 * (OutliersR.resCG Cq r νv * η) < ρv l := by
      intro l
      have h1 := hgb l
      linarith
    have hzlt : ∀ l : Fin r,
        lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) < ρv l := by
      intro l
      have h1 := hzedge l
      linarith
    have hE1 : ∀ (j l : Fin r), |R4.cform (m.W0hetR w N ω) (ρv l)
        (fun i => m.QmatHetR w N ω i j) (fun i => (m.QmatHetR w N ω * Z) i l)
        + Z j l| ≤ η := by
      intro j l
      have h := (hgood4 (j, l)).le
      rwa [sub_neg_eq_add] at h
    have hE2 : ∀ l : Fin r, |R4.cform2 (m.W0hetR w N ω) (ρv l)
        (fun i => (m.QmatHetR w N ω * Z) i l) (fun i => (m.QmatHetR w N ω * Z) i l)
        - νv l| ≤ η := fun l => (hgood5 l).le
    have hres := fun l => OutliersR.residual_bound_Z (m.isHermitian_W0hetR w N ω)
      (m.gram_eq_hetR w N ω) hνpos hzlt hCq0 hη0.le hcolb hE1 l
    have hgram := OutliersR.gram_bound_Z (m.isHermitian_W0hetR w N ω) hρinj hνpos hzlt
      hZb0 hη0.le hZbcol hE1 hE2
    have hbelow : ∀ l : Fin r, s ≤ (l : ℕ) →
        ρv l + 2 * (OutliersR.resCG Cq r νv * η) ≤ τ := by
      intro l hl
      have h1 := hbelow0 l hl
      linarith
    have hsep' : ∀ j l : Fin r, s ≤ (j : ℕ) → s ≤ (l : ℕ) → j ≠ l →
        4 * (OutliersR.resCG Cq r νv * η) < |ρv j - ρv l| := by
      intro j l _ _ hjl
      have h1 := hgsep j l hjl
      linarith
    have hI := Frame.count_eq_of_frame_of_split (m.isHermitian_W0hetR w N ω)
      (isHermitian_mul_transpose_self (m.stackXW w N ω)) (m.gram_eq_hetR w N ω) hsr
      hlamτ hmg hδ'0 hδ'half hzedge hbelow habove hsep' hres hgram hsmallc
    have hE1s : ∀ (j : Fin r) (a : Fin s), |R4.cform (m.W0hetR w N ω) (ρs a)
        (fun i => m.QmatHetR w N ω i j) (fun i => (m.QmatHetR w N ω * Zs) i a)
        + Zs j a| ≤ η := fun j a => hE1 j (Fin.castLE hsr a)
    have hE2s : ∀ a : Fin s, |R4.cform2 (m.W0hetR w N ω) (ρs a)
        (fun i => (m.QmatHetR w N ω * Zs) i a) (fun i => (m.QmatHetR w N ω * Zs) i a)
        - νs a| ≤ η := fun a => hE2 (Fin.castLE hsr a)
    have hE3s : ∀ a : Fin s, |R4.cform (m.W0hetR w N ω) (ρs a)
        ((Real.sqrt ((fun q => m.QmatHetR w N ω q k) ⬝ᵥ (fun q => m.QmatHetR w N ω q k)))⁻¹
          • fun q => m.QmatHetR w N ω q k)
        (fun i => (m.QmatHetR w N ω * Zs) i a) - tv a| ≤ η := by
      intro a
      rw [R4.cform_smul_left]
      exact (hgood6 a).le
    have hdet := OutliersR.align_detZ (m.isHermitian_W0hetR w N ω)
      (isHermitian_mul_transpose_self (m.stackXW w N ω)) (m.gram_eq_hetR w N ω) hρsinj hmg
      hlamτ hρsτ hνspos hI hCq0 hcolb hZb0 hZbs (OutliersR.dot_self_normalize hgpos) hLbt
      hη0 hη1 hE1s hE2s hE3s hsmalls
    have hdet' : |‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ) (Set.Ioi τ)
          (WithLp.toLp 2 ((Real.sqrt ((fun q => m.QmatHetR w N ω q k)
            ⬝ᵥ (fun q => m.QmatHetR w N ω q k)))⁻¹
              • fun q => m.QmatHetR w N ω q k))‖ ^ 2 - Tv| ≤ B * η :=
      le_trans hdet hRHS
    have hclose := OutliersR.abs_mul_sub_le_of_close hgpos hgle hdet' hgood3.le
    have hω' : ε ≤ |‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ) (Set.Ioi τ)
        (WithLp.toLp 2 fun q => m.QmatHetR w N ω q k)‖ ^ 2
        - ∑ l ∈ Finset.univ.filter (fun l : Fin r => (l : ℕ) < s),
            (WithLp.ofLp (z l) k) ^ 2
              / SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w (γ l) (z l)| := hω
    rw [OutliersR.normSq_specProj_eq_dot_mul _ _ hgpos, ← htarget] at hω'
    have hKe : (m.colGramLimitSW w c k + 1) * (B * η) = K * η := by rw [hKdef]; ring
    rw [hKe] at hclose
    linarith
  · exact tendsto_measure_zero_union hedgeC
      (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hcolT)
        (tendsto_measure_zero_union hcol2T
          (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE1T)
            (tendsto_measure_zero_union (tendsto_measure_zero_iUnion hE2T)
              (tendsto_measure_zero_iUnion hE3T)))))

end UnalignedModelR

end StackedSVD
