/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.Het.Align
import StackedSVD.RankR.Het.Bulk
import StackedSVD.RankR.Het.Simplicity

/-!
# F8, route D: drop the tables that do not carry component `j`

F8 (2026-09-05) weakened the positivity field of `SpikedModelR` to `hθnn : ∀ k, 0 ≤ θ k`. A
table may then have a zero spike at its last component. The stack weight of component `j` is
`w_ij = θ_ij/√(θ_ij² + c_i)`, which is `0` when `θ_ij = 0`. Table `i` then contributes
nothing to `X_stack^(j)`. The Gaussian chain of Track E needs a nonzero weight in every
table, so it does not apply to `m`.

This file is the drop step of route D (`notes/archive/F8_zero_spikes.md`, choice 4). It builds the
sub-model on the tables that carry component `j`, runs the per-component chain there, and
transports the conclusion back. The transport is exact, not approximate: the two weighted
Gram matrices are equal, because the dropped terms of `∑_i w_ij² X_iᵀ X_i` are zero.

## Contents

1. `UnalignedModelR.sub`: the sub-model along an embedding `e : Fin M' ↪ Fin M`, with the
   bridge lemmas, the regime, the alignment and the Gaussian marginal.
2. `stackGramW_sub`, `overlapIdx_congr_gram`, `SimpleIdx.congr_mat`: the transport tools.
3. `Scalars.gW_precomp` and `Scalars.stackSVDLimitW_precomp`: the scalar `γ_j` does not see
   a table whose strength at `j` is zero.
4. `key_of_gaussian_pos`: the per-component chain, moved here from `RankR/Het/Sup.lean`.
5. `key_of_gaussian` and `simpleIdxJ_of_gaussian_of_exists`: the same two conclusions under
   `hex : ∃ i, 0 < (m.tbl i).θ j` alone.

No `sorry`, no `axiom`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

/-! ### 1. Transport tools for one index -/

/-- Two matrices with the same Gram matrix have the same index overlaps. The row counts may
differ, which is what the drop step needs. -/
theorem overlapIdx_congr_gram {p q q' : ℕ} (X : Matrix (Fin q) (Fin p) ℝ)
    (X' : Matrix (Fin q') (Fin p) ℝ) (h : Xᵀ * X = X'ᵀ * X') (k : ℕ)
    (w : EuclideanSpace ℝ (Fin p)) : overlapIdx X k w = overlapIdx X' k w := by
  unfold overlapIdx
  rw [specProjIdx_congr_mat (isHermitian_transpose_mul_self X)
    (isHermitian_transpose_mul_self X') h k]

/-- `SimpleIdx` only reads the matrix, so equal matrices give equivalent statements. -/
theorem SimpleIdx.congr_mat {p : ℕ} {A B : Matrix (Fin p) (Fin p) ℝ} (hA : A.IsHermitian)
    (hB : B.IsHermitian) (h : A = B) (k : ℕ) : SimpleIdx A hA k ↔ SimpleIdx B hB k := by
  subst h; exact Iff.rfl

/-- The image of an embedding, as a `Finset` membership. -/
theorem mem_map_univ_of_mem_range {α β : Type*} [Fintype α] {e : α ↪ β} {i : β}
    (h : i ∈ Set.range e) : i ∈ Finset.univ.map e := by
  obtain ⟨a, rfl⟩ := h
  exact Finset.mem_map_of_mem e (Finset.mem_univ a)

/-! ### 2. The scalars do not see a zero strength -/

namespace Scalars

variable {M M' r : ℕ}

/-- A table with strength `0` at component `j` gets stack weight `0`. -/
theorem wStackR_eq_zero {θ : Fin M → Fin r → ℝ} {c : Fin M → ℝ} {j : Fin r} {i : Fin M}
    (h : θ i j = 0) : wStackR θ c j i = 0 := by
  rw [wStackR_apply, h, zero_div]

/-- The secular function of `stackSVDLimitW` drops the tables with strength `0`. -/
theorem gW_precomp (θ c : Fin M → ℝ) (e : Fin M' ↪ Fin M)
    (hz : ∀ i, i ∉ Set.range e → θ i = 0) : gW (θ ∘ e) (c ∘ e) = gW θ c := by
  funext x
  rw [gW, gW]
  simp only [Function.comp_apply]
  rw [← Finset.sum_map Finset.univ e
    (fun i => θ i ^ 4 * (1 - x) / (c i + x * θ i ^ 2))]
  refine Finset.sum_subset (Finset.subset_univ _) fun i _ hi => ?_
  rw [hz i fun hr => hi (mem_map_univ_of_mem_range hr)]
  norm_num

/-- `γ` does not see a table whose strength is `0`: the root of the secular equation is the
same on the sub-family. -/
theorem stackSVDLimitW_precomp (θ c : Fin M → ℝ) (e : Fin M' ↪ Fin M)
    (hz : ∀ i, i ∉ Set.range e → θ i = 0) :
    stackSVDLimitW (θ ∘ e) (c ∘ e) = stackSVDLimitW θ c := by
  have hg : gW (θ ∘ e) (c ∘ e) = gW θ c := gW_precomp θ c e hz
  unfold stackSVDLimitW
  by_cases hex : ∃! x : ℝ, x ∈ Set.Ioo (0:ℝ) 1 ∧ gW θ c x = 1
  · have hex' : ∃! x : ℝ, x ∈ Set.Ioo (0:ℝ) 1 ∧ gW (θ ∘ e) (c ∘ e) x = 1 := by
      rw [hg]; exact hex
    rw [dif_pos hex, dif_pos hex']
    refine hex.choose_spec.2 _ ⟨hex'.choose_spec.1.1, ?_⟩
    rw [← hg]; exact hex'.choose_spec.1.2
  · have hex' : ¬ ∃! x : ℝ, x ∈ Set.Ioo (0:ℝ) 1 ∧ gW (θ ∘ e) (c ∘ e) x = 1 := by
      rw [hg]; exact hex
    rw [dif_neg hex, dif_neg hex']

end Scalars

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ}

/-! ### 3. The sub-model -/

/-- The sub-model on the tables `e 0, …, e (M' - 1)`: the same shared `V`, with the tables
and the alignment matrices pulled back along `e`. -/
def sub (m : UnalignedModelR μ M n d r (alignedRk M r)) {M' : ℕ} (e : Fin M' ↪ Fin M) :
    UnalignedModelR μ M' (fun i' => n (e i')) d r (alignedRk M' r) where
  tbl := fun i' => m.tbl (e i')
  V := m.V
  R := fun i' => m.R (e i')
  hV := m.hV
  hR := fun i' => m.hR (e i')
  hv := fun i' N => m.hv (e i') N

section Sub

variable (m : UnalignedModelR μ M n d r (alignedRk M r)) {M' : ℕ} (e : Fin M' ↪ Fin M)

@[simp] theorem sub_tbl (i' : Fin M') : (m.sub e).tbl i' = m.tbl (e i') := rfl

@[simp] theorem sub_V (N : ℕ) : (m.sub e).V N = m.V N := rfl

@[simp] theorem sub_R (i' : Fin M') : (m.sub e).R i' = m.R (e i') := rfl

@[simp] theorem sub_colVecG (N : ℕ) (k : Fin r) : (m.sub e).colVecG N k = m.colVecG N k := rfl

@[simp] theorem sub_thetaAligned :
    (m.sub e).thetaAligned = fun i' k => m.thetaAligned (e i') k := rfl

/-- The regime of every kept table is inherited. -/
theorem sub_Regime {c : Fin M → ℝ} (hreg : ∀ i, (m.tbl i).Regime (c i)) :
    ∀ i', ((m.sub e).tbl i').Regime (c (e i')) := fun i' => hreg (e i')

/-- Exact alignment is inherited. -/
theorem sub_hR (hR : ∀ i, m.R i = 1) : ∀ i', (m.sub e).R i' = 1 := fun i' => hR (e i')

/-- **The Gaussian marginal.** The joint law of the kept noise tables is the product of their
laws: the full family is independent with Gaussian marginals, and `iIndepFun.precomp` keeps
both properties along an injection. -/
theorem sub_JointGaussianNoise (hG : m.JointGaussianNoise) : (m.sub e).JointGaussianNoise := by
  intro N
  have : IsProbabilityMeasure (μ N) := (hG N).isProbabilityMeasure
  have h1 : ∀ i, HasLaw (fun ω => (m.tbl i).Z N ω) (gaussianMatrix (n i N) (d N)) (μ N) :=
    fun i => (MeasureTheory.measurePreserving_eval
      (fun i : Fin M => gaussianMatrix (n i N) (d N)) i).fun_comp_hasLaw (hG N)
  have h2 : iIndepFun (fun i ω => (m.tbl i).Z N ω) (μ N) :=
    (iIndepFun_iff_hasLaw_pi_pi h1).mpr (hG N)
  exact (h2.precomp e.injective).hasLaw_pi fun i' => h1 (e i')

/-- **The Gram identity.** With `w = 0` off the range of `e`, the weighted Gram matrix of `m`
is the weighted Gram matrix of the sub-model. -/
theorem stackGramW_sub (w : Fin M → ℝ) (hw : ∀ i, i ∉ Set.range e → w i = 0) (N : ℕ)
    (ω : Ω N) : m.stackGramW w N ω = (m.sub e).stackGramW (fun i' => w (e i')) N ω := by
  rw [stackGramW_eq_sum, stackGramW_eq_sum]
  have hmap : ∑ i ∈ Finset.univ.map e,
        w i ^ 2 • (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω)
      = ∑ i, w i ^ 2 • (((m.tbl i).X N ω)ᵀ * (m.tbl i).X N ω) := by
    refine Finset.sum_subset (Finset.subset_univ _) fun i _ hi => ?_
    rw [hw i fun hr => hi (mem_map_univ_of_mem_range hr)]
    simp
  rw [← hmap, Finset.sum_map]
  rfl

/-- The stack weights of the sub-model are the stack weights of `m` at the kept tables. -/
theorem wStackR_sub (c : Fin M → ℝ) (j : Fin r) :
    Scalars.wStackR (m.sub e).thetaAligned (fun i' => c (e i')) j
      = fun i' => Scalars.wStackR m.thetaAligned c j (e i') := rfl

/-- `γ_j` is the same on the sub-model when every dropped table has strength `0` at `j`. -/
theorem gammaR_sub (c : Fin M → ℝ) (j : Fin r)
    (hz : ∀ i, i ∉ Set.range e → (m.tbl i).θ j = 0) :
    Scalars.gammaR (m.sub e).thetaAligned (fun i' => c (e i')) j
      = Scalars.gammaR m.thetaAligned c j :=
  Scalars.stackSVDLimitW_precomp (fun i => m.thetaAligned i j) c e hz

end Sub

/-! ### 4. The per-component chain, with every table carrying component `j` -/

/-- **The two limit fields of `HeteroLawR` at one component, under `hposj`.** Read the
squared overlap of column `l` of the signal frame with the `ℓ_j`-th right singular subspace
of `X_stack^(j)`. It tends to `γ_j` when `l = j`, and to `0` when `l ≠ j`.

The proof splits on `eq:assumption4` at component `j`, at the weights of component `j`.
Above the threshold stage E8a (`align_cross_het_of_gaussian_sup`) applies at the index
`ellSup`, which the class identifies with `ℓ_j` (`Scalars.ellSup_eq_ellR`). Below it
`γ_j = 0` (`Scalars.gammaR_eq_zero_of_not_assumption4`) and `ℓ_j` is at or above `numSup`
(`Scalars.numSup_le_of_not_assumption4`), so stage E8b (`cross_het_of_gaussian_bulk`) gives
the limit `0` at every column. This block was the `key` of `heteroLawR_of_gaussian_aux`
before F8. -/
theorem key_of_gaussian_pos [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r) {j : Fin r} (hposj : ∀ i, 0 < (m.tbl i).θ j)
    (l : Fin r) :
    TendstoInProb μ
      (fun N ω => overlapIdx (m.stackXJ c j N ω) (Scalars.ellR m.thetaAligned c j)
        (m.colVecG N l))
      (if l = j then Scalars.gammaR m.thetaAligned c j else 0) := by
  have hnn : ∀ i k, 0 ≤ m.thetaAligned i k := fun i k => (m.tbl i).hθnn k
  have hanti : ∀ i, StrictAnti (m.thetaAligned i) := fun i => (m.tbl i).hθanti
  have hM : 0 < M := NeZero.pos M
  have hw : ∀ i, Scalars.wStackR m.thetaAligned c j i ≠ 0 :=
    fun i => m.wStackR_thetaAligned_ne_zero hc hposj i
  have hedge := m.heteroEdgeR_of_gaussian (Scalars.wStackR m.thetaAligned c j) c hc
    ⟨⟨0, hM⟩, hw _⟩ hreg hG hpd
  by_cases h4 : Scalars.Assumption4 (fun i => m.thetaAligned i j) c
      (Scalars.wStackR m.thetaAligned c j)
  · -- Supercritical at component `j`: stage E8a at the index `ellSup … j = ℓ_j`.
    have h := m.align_cross_het_of_gaussian_sup (Scalars.wStackR m.thetaAligned c j) c hc
      hw hR hreg hG hpd hedge h4 l
    rw [Scalars.ellSup_eq_ellR hc hnn hanti hM hposj h4,
      Scalars.Lw_wStackR_eq_gammaR hc j] at h
    exact h
  · -- Subcritical at component `j`: `γ_j = 0` and `ℓ_j` is a bulk index, so stage E8b.
    rw [Scalars.gammaR_eq_zero_of_not_assumption4 hc h4, ite_self]
    exact m.cross_het_of_gaussian_bulk (Scalars.wStackR m.thetaAligned c j) c hc hw hR hreg
      hG hpd hedge (a := ⟨Scalars.ellR m.thetaAligned c j,
        Scalars.ellR_lt m.thetaAligned c j⟩)
      (Scalars.numSup_le_of_not_assumption4 hc hnn hanti hM hposj h4) l

/-! ### 5. The drop step -/

section Carrier

variable (m : UnalignedModelR μ M n d r (alignedRk M r)) (j : Fin r)

/-- The tables that carry component `j`, that is, the tables with `θ_ij > 0`. -/
noncomputable def carrier : Finset (Fin M) := Finset.univ.filter fun i => 0 < (m.tbl i).θ j

/-- The carrier of component `j`, as an embedding out of `Fin (carrier).card`. -/
noncomputable def carrierEmb : Fin (m.carrier j).card ↪ Fin M :=
  ((m.carrier j).orderEmbOfFin rfl).toEmbedding

/-- The range of `carrierEmb` is exactly the set of tables that carry component `j`. -/
theorem mem_range_carrierEmb (i : Fin M) :
    i ∈ Set.range (m.carrierEmb j) ↔ 0 < (m.tbl i).θ j := by
  change i ∈ Set.range ((m.carrier j).orderEmbOfFin rfl) ↔ _
  rw [Finset.range_orderEmbOfFin]
  simp [carrier]

/-- A dropped table has strength `0` at component `j`, by `hθnn`. -/
theorem theta_eq_zero_of_not_mem_range {i : Fin M} (hi : i ∉ Set.range (m.carrierEmb j)) :
    (m.tbl i).θ j = 0 :=
  le_antisymm (not_lt.mp fun h => hi ((m.mem_range_carrierEmb j i).mpr h)) ((m.tbl i).hθnn j)

/-- A dropped table has stack weight `0` at component `j`. -/
theorem wStackR_eq_zero_of_not_mem_range (c : Fin M → ℝ) {i : Fin M}
    (hi : i ∉ Set.range (m.carrierEmb j)) : Scalars.wStackR m.thetaAligned c j i = 0 :=
  Scalars.wStackR_eq_zero (m.theta_eq_zero_of_not_mem_range j hi)

/-- Every kept table carries component `j`. -/
theorem carrier_pos (i' : Fin (m.carrier j).card) :
    0 < ((m.sub (m.carrierEmb j)).tbl i').θ j :=
  (m.mem_range_carrierEmb j _).mp ⟨i', rfl⟩

/-- One table with `θ_ij > 0` makes the carrier nonempty. -/
theorem neZero_carrier_card (hex : ∃ i, 0 < (m.tbl i).θ j) : NeZero (m.carrier j).card := by
  obtain ⟨i₀, hi₀⟩ := hex
  exact ⟨(Finset.card_pos.mpr ⟨i₀, by simp [carrier, hi₀]⟩).ne'⟩

/-- **The Gram matrices agree.** The dropped tables have weight `0`, so they contribute
nothing to `∑_i w_ij² X_iᵀ X_i`. -/
theorem stackGramJ_sub_carrier (c : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    m.stackGramJ c j N ω
      = (m.sub (m.carrierEmb j)).stackGramJ (fun i' => c (m.carrierEmb j i')) j N ω := by
  change m.stackGramW (Scalars.wStackR m.thetaAligned c j) N ω
    = (m.sub (m.carrierEmb j)).stackGramW
        (Scalars.wStackR (m.sub (m.carrierEmb j)).thetaAligned
          (fun i' => c (m.carrierEmb j i')) j) N ω
  rw [m.wStackR_sub (m.carrierEmb j) c j]
  exact m.stackGramW_sub (m.carrierEmb j) _
    (fun i hi => m.wStackR_eq_zero_of_not_mem_range j c hi) N ω

/-- Both models read the same index: it is `j` on each side (`ellR_thetaAligned`). -/
theorem ellR_sub_carrier {c : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hex : ∃ i, 0 < (m.tbl i).θ j) :
    Scalars.ellR (m.sub (m.carrierEmb j)).thetaAligned (fun i' => c (m.carrierEmb j i')) j
      = Scalars.ellR m.thetaAligned c j := by
  have := m.neZero_carrier_card j hex
  rw [(m.sub (m.carrierEmb j)).ellR_thetaAligned _ (fun i' => hc _) ⟨0, m.carrier_pos j 0⟩,
    m.ellR_thetaAligned c hc hex]

end Carrier

/-- **Route D: the component `j` limits under `∃ i, 0 < θ_ij` alone.** The tables with
`θ_ij = 0` carry weight `0` in `X_stack^(j)`. The proof drops them through the sub-model on
the carrier of `j`, where `key_of_gaussian_pos` applies. The transport is exact: the two
weighted Gram matrices are equal, `ℓ_j` is `j` on both sides, and `γ_j` does not see a table
of strength `0` (`Scalars.stackSVDLimitW_precomp`). -/
theorem key_of_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r) {j : Fin r} (hex : ∃ i, 0 < (m.tbl i).θ j)
    (l : Fin r) :
    TendstoInProb μ
      (fun N ω => overlapIdx (m.stackXJ c j N ω) (Scalars.ellR m.thetaAligned c j)
        (m.colVecG N l))
      (if l = j then Scalars.gammaR m.thetaAligned c j else 0) := by
  have := m.neZero_carrier_card j hex
  set e := m.carrierEmb j with he
  have h := (m.sub e).key_of_gaussian_pos (fun i' => c (e i')) (fun i' => hc _)
    (m.sub_hR e hR) (m.sub_Regime e hreg) (m.sub_JointGaussianNoise e hG) hpd
    (m.carrier_pos j) l
  rw [m.ellR_sub_carrier j hc hex,
    m.gammaR_sub e c j (fun i hi => m.theta_eq_zero_of_not_mem_range j hi)] at h
  have hfun : (fun N (ω : Ω N) => overlapIdx (m.stackXJ c j N ω)
        (Scalars.ellR m.thetaAligned c j) (m.colVecG N l))
      = fun N ω => overlapIdx ((m.sub e).stackXJ (fun i' => c (e i')) j N ω)
        (Scalars.ellR m.thetaAligned c j) ((m.sub e).colVecG N l) := by
    funext N ω
    exact overlapIdx_congr_gram _ _ (m.stackGramJ_sub_carrier j c N ω) _ _
  rw [hfun]
  exact h

/-- **`simpleIdxJ` under `∃ i, 0 < θ_ij` alone.** The same drop step: the Gram matrix of
`X_stack^(j)` is the Gram matrix of the sub-model on the carrier of `j`, and `SimpleIdx`
reads only that matrix. -/
theorem simpleIdxJ_of_gaussian_of_exists [NeZero M]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) {c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hG : m.JointGaussianNoise) {j : Fin r} (hex : ∃ i, 0 < (m.tbl i).θ j) (N : ℕ) :
    ∀ᵐ ω ∂(μ N), SimpleIdx (m.stackGramJ c j N ω) (m.isHermitian_stackGramJ c j N ω)
      (Scalars.ellR m.thetaAligned c j) := by
  have := m.neZero_carrier_card j hex
  set e := m.carrierEmb j with he
  have h := (m.sub e).simpleIdxJ_of_gaussian (c := fun i' => c (e i')) (fun i' => hc _)
    (m.sub_JointGaussianNoise e hG) (m.carrier_pos j) N
  rw [m.ellR_sub_carrier j hc hex] at h
  filter_upwards [h] with ω hω
  exact (SimpleIdx.congr_mat _ _ (m.stackGramJ_sub_carrier j c N ω) _).mpr hω

end UnalignedModelR

end StackedSVD
