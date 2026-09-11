/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.Het.Sub
import StackedSVD.RankR.Het.Simplicity
import StackedSVD.RankR.SingleWeight.Scalars

/-!
# Track G, item F18c: drop the tables of weight zero

`prop:gen_rank_stacksvd_singleweight` (`main_paper.tex:2112`) puts no condition on the
weights `w_i`. The Gaussian chain of Track G met the paper's statement except in one place:
the almost sure simplicity of the weighted stack Gram matrix (`simpleSpec_ae_stackGramW`,
`RankR/Het/Simplicity.lean`) needs a nonzero weight in every table, and the first version of
`RankR/SingleWeight/Het/Sup.lean` (2026-09-05) took `hw : ∀ i, w i ≠ 0` for that reason
(decision D37 of `notes/FLAGGED.md`).

This file is the drop step that removes `hw`, the twin of route D of F8
(`RankR/Het/Sub.lean`) at general `rk`. A table of weight zero contributes nothing to
`X_stack(w)ᵀ X_stack(w) = ∑_i w_i² X_iᵀ X_i`, so the Gram matrix of `m` is the Gram matrix of
the sub-model on the tables of nonzero weight, where every weight is nonzero. `SimpleSpec`
reads only that matrix. The transport is exact.

## Contents

1. `SingleWeight.EigSep.exists_ne_zero`: `assum:gen_rank_stacksvd_eig_sep` itself forces a
   nonzero weight once `r ≥ 1`, since the secular matrix of the zero weight is `0` and
   `det (1 - 0) = 1`. So no hypothesis on `w` remains in the Gaussian theorems.
2. `UnalignedModelR.subRk`: the sub-model along an embedding `e : Fin M' ↪ Fin M` at general
   `rk` (`UnalignedModelR.sub` of `RankR/Het/Sub.lean` is fixed at `alignedRk M r`), with
   the regime, the Gaussian marginal and the Gram identity `stackGramW_subRk`.
3. `wSupp`, `wSuppEmb`: the support of `w` as a `Finset` and as an embedding, and
   `sum_wSuppEmb`, the sum over the support.
4. `simpleSpec_ae_stackGramW_of_exists`: simplicity under `∃ i, w i ≠ 0`, with the size
   condition read on the support, `k ≤ min (∑ i with w i ≠ 0, n i N) (d N)`.

No `sorry`, no `axiom`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

/-! ### 1. `EigSep` forces a nonzero weight -/

namespace SingleWeight

variable {M r : ℕ} {rk : Fin M → ℕ} {θ : (i : Fin M) → Fin (rk i) → ℝ}
  {R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ} {w c : Fin M → ℝ}

/-- The secular matrix of the zero weight is `0`. -/
theorem secMat_eq_zero_of_forall (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (hw : ∀ i, w i = 0) (γ : ℝ) :
    secMat θ R w γ = 0 := by
  unfold secMat
  simp [hw]

/-- **`EigSep` needs a nonzero weight.** At the zero weight the secular equation reads
`det (1 - 0) = 0`, which is false, so `IsSecularRoot` fails at every `γ`. -/
theorem EigSep.exists_ne_zero {γ : Fin r → ℝ} {z : Fin r → EuclideanSpace ℝ (Fin r)}
    (hsep : EigSep θ R w c γ z) (l : Fin r) : ∃ i, w i ≠ 0 := by
  by_contra h
  push Not at h
  have h1 := (hsep.root l).2
  rw [secMat_eq_zero_of_forall θ R h, sub_zero, Matrix.det_one] at h1
  exact one_ne_zero h1

/-- Either `r = 0` (every statement about the components is vacuous) or `EigSep` supplies a
table of nonzero weight. -/
theorem EigSep.eq_zero_or_exists {γ : Fin r → ℝ} {z : Fin r → EuclideanSpace ℝ (Fin r)}
    (hsep : EigSep θ R w c γ z) : r = 0 ∨ ∃ i, w i ≠ 0 := by
  rcases Nat.eq_zero_or_pos r with hr | hr
  · exact Or.inl hr
  · exact Or.inr (hsep.exists_ne_zero ⟨0, hr⟩)

end SingleWeight

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-! ### 2. The sub-model at general `rk` -/

/-- The sub-model on the tables `e 0, …, e (M' - 1)` at general `rk`: the same shared `V`,
with the tables, the ranks and the alignment matrices pulled back along `e`. Twin of
`UnalignedModelR.sub` (`RankR/Het/Sub.lean:116`), which is fixed at `alignedRk M r`. -/
def subRk (m : UnalignedModelR μ M n d r rk) {M' : ℕ} (e : Fin M' ↪ Fin M) :
    UnalignedModelR μ M' (fun i' => n (e i')) d r (fun i' => rk (e i')) where
  tbl := fun i' => m.tbl (e i')
  V := m.V
  R := fun i' => m.R (e i')
  hV := m.hV
  hR := fun i' => m.hR (e i')
  hv := fun i' N => m.hv (e i') N

section SubRk

variable (m : UnalignedModelR μ M n d r rk) {M' : ℕ} (e : Fin M' ↪ Fin M)

@[simp] theorem subRk_tbl (i' : Fin M') : (m.subRk e).tbl i' = m.tbl (e i') := rfl

@[simp] theorem subRk_V (N : ℕ) : (m.subRk e).V N = m.V N := rfl

@[simp] theorem subRk_R (i' : Fin M') : (m.subRk e).R i' = m.R (e i') := rfl

/-- The regime of every kept table is inherited. -/
theorem subRk_Regime {c : Fin M → ℝ} (hreg : ∀ i, (m.tbl i).Regime (c i)) :
    ∀ i', ((m.subRk e).tbl i').Regime (c (e i')) := fun i' => hreg (e i')

/-- **The Gaussian marginal.** Copy of `sub_JointGaussianNoise` (`RankR/Het/Sub.lean:150`):
the full family is independent with Gaussian marginals, and `iIndepFun.precomp` keeps both
properties along an injection. -/
theorem subRk_JointGaussianNoise (hG : m.JointGaussianNoise) :
    (m.subRk e).JointGaussianNoise := by
  intro N
  have : IsProbabilityMeasure (μ N) := (hG N).isProbabilityMeasure
  have h1 : ∀ i, HasLaw (fun ω => (m.tbl i).Z N ω) (gaussianMatrix (n i N) (d N)) (μ N) :=
    fun i => (MeasureTheory.measurePreserving_eval
      (fun i : Fin M => gaussianMatrix (n i N) (d N)) i).fun_comp_hasLaw (hG N)
  have h2 : iIndepFun (fun i ω => (m.tbl i).Z N ω) (μ N) :=
    (iIndepFun_iff_hasLaw_pi_pi h1).mpr (hG N)
  exact (h2.precomp e.injective).hasLaw_pi fun i' => h1 (e i')

/-- **The Gram identity.** With `w = 0` off the range of `e`, the weighted Gram matrix of `m`
is the weighted Gram matrix of the sub-model. Copy of `stackGramW_sub`
(`RankR/Het/Sub.lean:162`). -/
theorem stackGramW_subRk (w : Fin M → ℝ) (hw : ∀ i, i ∉ Set.range e → w i = 0) (N : ℕ)
    (ω : Ω N) : m.stackGramW w N ω = (m.subRk e).stackGramW (fun i' => w (e i')) N ω := by
  rw [stackGramW_eq_sum, stackGramW_eq_sum]
  have hmap : ∑ i ∈ Finset.univ.map e,
        w i ^ 2 • (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω)
      = ∑ i, w i ^ 2 • (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω) := by
    refine Finset.sum_subset (Finset.subset_univ _) fun i _ hi => ?_
    rw [hw i fun hr => hi (mem_map_univ_of_mem_range hr)]
    simp
  rw [← hmap, Finset.sum_map]
  rfl

end SubRk

/-! ### 3. The support of the weights -/

section Supp

variable (w : Fin M → ℝ)

/-- The tables of nonzero weight. -/
noncomputable def wSupp : Finset (Fin M) := Finset.univ.filter fun i => w i ≠ 0

/-- The support of `w`, as an embedding out of `Fin (wSupp w).card`. -/
noncomputable def wSuppEmb : Fin (wSupp w).card ↪ Fin M :=
  ((wSupp w).orderEmbOfFin rfl).toEmbedding

/-- The range of `wSuppEmb` is exactly the set of tables of nonzero weight. -/
theorem mem_range_wSuppEmb (i : Fin M) : i ∈ Set.range (wSuppEmb w) ↔ w i ≠ 0 := by
  change i ∈ Set.range ((wSupp w).orderEmbOfFin rfl) ↔ _
  rw [Finset.range_orderEmbOfFin]
  simp [wSupp]

/-- A dropped table has weight `0`. -/
theorem eq_zero_of_not_mem_range_wSuppEmb {i : Fin M} (hi : i ∉ Set.range (wSuppEmb w)) :
    w i = 0 :=
  not_not.mp fun h => hi ((mem_range_wSuppEmb w i).mpr h)

/-- Every kept table has nonzero weight. -/
theorem wSuppEmb_ne_zero (i' : Fin (wSupp w).card) : w (wSuppEmb w i') ≠ 0 :=
  (mem_range_wSuppEmb w _).mp ⟨i', rfl⟩

/-- One table of nonzero weight makes the support nonempty. -/
theorem neZero_wSupp_card (hex : ∃ i, w i ≠ 0) : NeZero (wSupp w).card := by
  obtain ⟨i₀, hi₀⟩ := hex
  exact ⟨(Finset.card_pos.mpr ⟨i₀, by simp [wSupp, hi₀]⟩).ne'⟩

/-- A sum over the support, read through the embedding. -/
theorem sum_wSuppEmb (f : Fin M → ℕ) :
    ∑ i', f (wSuppEmb w i') = ∑ i with w i ≠ 0, f i := by
  rw [← Finset.sum_map, wSuppEmb, Finset.map_orderEmbOfFin_univ]
  rfl

end Supp

/-! ### 4. Simplicity under `∃ i, w i ≠ 0` -/

/-- **The weighted stack Gram matrix has pairwise distinct top-`k` eigenvalues almost surely,
when some weight is nonzero.** The size condition reads the rows of the tables of nonzero
weight, since the tables of weight zero contribute nothing to the Gram matrix. Proof: the
Gram matrix of `m` is the Gram matrix of the sub-model on the support
(`stackGramW_subRk`), where `simpleSpec_ae_stackGramW` (`RankR/Het/Simplicity.lean:161`)
applies. -/
theorem simpleSpec_ae_stackGramW_of_exists (m : UnalignedModelR μ M n d r rk)
    (w : Fin M → ℝ) (hG : m.JointGaussianNoise) (hex : ∃ i, w i ≠ 0) (k N : ℕ)
    (hk : k ≤ min (∑ i with w i ≠ 0, n i N) (d N)) :
    ∀ᵐ ω ∂(μ N), SimpleSpec (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) k := by
  have := neZero_wSupp_card w hex
  set e := wSuppEmb w with he
  have hsum : ∑ i', n (e i') N = ∑ i with w i ≠ 0, n i N := sum_wSuppEmb w fun i => n i N
  have hk' : k ≤ min (∑ i', n (e i') N) (d N) := by rw [hsum]; exact hk
  have h := (m.subRk e).simpleSpec_ae_stackGramW (fun i' => w (e i'))
    (m.subRk_JointGaussianNoise e hG) (fun i' => wSuppEmb_ne_zero w i') k N hk'
  filter_upwards [h] with ω hω
  exact EdgeGlueR.simpleSpec_congr
    (m.stackGramW_subRk e w (fun i hi => eq_zero_of_not_mem_range_wSuppEmb w hi) N ω).symm
    _ _ hω

end UnalignedModelR

end StackedSVD
