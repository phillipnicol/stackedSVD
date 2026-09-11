/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.SingleWeight.Het.Scalars
import StackedSVD.RankR.SingleWeight.Het.Outliers
import StackedSVD.RankR.Het.Align
import StackedSVD.RankR.Het.Duality
import StackedSVD.RankR.Het.Deloc
import StackedSVD.LinAlg.SpecWindow

/-!
# Track G, unit G4: the eigenvalue limit and the `align` field of `SingleWeightLaw`

Unit G4 of `notes/archive/trackG_plan.md` (2026-09-05). Mirrors `RankR/Het/Align.lean:322`
(`tendstoInProb_eigVal_het`) and `:444` (`align_cross_het_of_gaussian_sup`) at general `R_i`,
under `SingleWeight.EigSep` instead of `Scalars.Assumption4`. The target
`align_sw_of_gaussian` is the `align` field of `SingleWeightLaw`
(`RankR/SingleWeight/Main.lean:67`).

## The route in five lines

1. Under `EigSep` every component is supercritical and `swRho c w (γ l)` is strictly
   decreasing in `l` (`EigSep.strictAnti_swRho`, unit G3), so the sorted index of outlier `l`
   is `l` itself. The aligned file's `Scalars.ellSup` counter disappears.
2. `exists_margin_sw` turns the margin `g` of unit G3 into a `δ > 0` that separates
   `ρ_l = swRho c w (γ l)` from the bulk edge and from every other outlier by `2 δ`.
3. At the thresholds `τ∓ = ρ_l ∓ δ` the gap conditions `gap_sub_of_margin` and
   `gap_add_of_margin` feed the count theorem of unit G3 at `s = l + 1` and `s = l`, which
   traps the sorted eigenvalue `λ_l` in `(τ₋, τ₊]` with probability tending to 1.
4. On that event the window lemmas `specProjIdx_eq_specProj_Ioc` and `normSq_specProj_Ioc`
   write the index projector norm as the difference of two half-line norms, and
   `sum_filter_lt_succ_sub` reduces the difference of the two half-line limits of unit G3 to
   the single term `(z_l)_k² / ν_l`.
5. The quotient by `λ_l → ρ_l` and the overlap identity `1/(ρ_l ν_l) = swTerm` (unit G0)
   give the limit `swTerm … (γ l) (z l) * (z_l)_k²` of the paper display at
   `main_paper.tex:2153`.

The index range `(l : ℕ) < ∑ i, n i N` is not a model constraint at general `rk`, so it
enters through `UnalignedModelR.eventually_r_le_min` and the eventual squeeze
`OutliersR.tendsto_measure_zero_of_eventually_subset`.

Every declaration is proved in full.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

/-! ### 1. Deterministic helpers -/

namespace SingleWeight

/-- **The margin.** From the margin `g` of `exists_margin_of_strictAnti` (unit G3) at
`δ = g / 4`: the family `ρ` sits `2 δ` above the floor `b`, and distinct members are `2 δ`
apart. The aligned mirror is `UnalignedModelR.exists_margin_het`
(`RankR/Het/Align.lean:76`), where the separation is a hypothesis; here the strict order of
the roots (`EigSep.sorted`) supplies it. -/
theorem exists_margin_sw {r : ℕ} {ρ : Fin r → ℝ} {b : ℝ}
    (hb : ∀ l, b < ρ l) (hinj : Function.Injective ρ) :
    ∃ δ > 0, (∀ l, b + 2 * δ < ρ l) ∧ ∀ k l, k ≠ l → 2 * δ < |ρ k - ρ l| := by
  obtain ⟨g, hg0, hgb, hgsep⟩ := exists_margin_of_strictAnti hb hinj
  refine ⟨g / 4, by linarith, fun l => ?_, fun k l hkl => ?_⟩
  · have := hgb l
    linarith
  · have := hgsep k l hkl
    linarith

/-- The gap hypothesis of unit G3 at the lower threshold `τ = ρ_l - δ` and margin `mg = δ / 2`
with `s = l + 1`: the outliers of index at most `l` sit at `τ + δ / 2` or above, and every
later outlier sits `δ / 2` below `τ`. -/
theorem gap_sub_of_margin {r : ℕ} {ρ : Fin r → ℝ} (hanti : StrictAnti ρ) {δ : ℝ}
    (hδ : 0 < δ) (hgap2 : ∀ k l, k ≠ l → 2 * δ < |ρ k - ρ l|) (l : Fin r) :
    ∀ l' : Fin r, ((l' : ℕ) < (l : ℕ) + 1 → (ρ l - δ) + δ / 2 ≤ ρ l')
      ∧ ((l : ℕ) + 1 ≤ (l' : ℕ) → ρ l' + δ / 2 < ρ l - δ) := by
  intro l'
  refine ⟨fun h => ?_, fun h => ?_⟩
  · rcases lt_or_ge (l' : ℕ) (l : ℕ) with hlt | hge
    · have := hanti (Fin.lt_def.mpr hlt)
      linarith
    · have hll : l' = l := Fin.val_injective (by omega)
      subst hll
      linarith
  · have hlt : l < l' := Fin.lt_def.mpr (by omega)
    have hlo := hanti hlt
    have hne : l ≠ l' := ne_of_lt hlt
    have h2 := hgap2 l l' hne
    rw [abs_of_pos (by linarith)] at h2
    linarith

/-- The gap hypothesis of unit G3 at the upper threshold `τ = ρ_l + δ` and margin
`mg = δ / 2` with `s = l`: the outliers of index below `l` sit at `τ + δ / 2` or above, and
every outlier of index at least `l` sits `δ / 2` below `τ`. -/
theorem gap_add_of_margin {r : ℕ} {ρ : Fin r → ℝ} (hanti : StrictAnti ρ) {δ : ℝ}
    (hδ : 0 < δ) (hgap2 : ∀ k l, k ≠ l → 2 * δ < |ρ k - ρ l|) (l : Fin r) :
    ∀ l' : Fin r, ((l' : ℕ) < (l : ℕ) → (ρ l + δ) + δ / 2 ≤ ρ l')
      ∧ ((l : ℕ) ≤ (l' : ℕ) → ρ l' + δ / 2 < ρ l + δ) := by
  intro l'
  refine ⟨fun h => ?_, fun h => ?_⟩
  · have hlt : l' < l := Fin.lt_def.mpr h
    have hlo := hanti hlt
    have hne : l' ≠ l := ne_of_lt hlt
    have h2 := hgap2 l' l hne
    rw [abs_of_pos (by linarith)] at h2
    linarith
  · rcases lt_or_ge (l : ℕ) (l' : ℕ) with hlt | hge
    · have := hanti (Fin.lt_def.mpr hlt)
      linarith
    · have hll : l' = l := Fin.val_injective (by omega)
      subst hll
      linarith

/-- The two half-line sums of unit G3 at `s = l + 1` and `s = l` differ by the single term
`f l`. -/
theorem sum_filter_lt_succ_sub {r : ℕ} (f : Fin r → ℝ) (l : Fin r) :
    (∑ l' ∈ Finset.univ.filter (fun l' : Fin r => (l' : ℕ) < (l : ℕ) + 1), f l')
      - ∑ l' ∈ Finset.univ.filter (fun l' : Fin r => (l' : ℕ) < (l : ℕ)), f l' = f l := by
  have hins : (Finset.univ.filter (fun l' : Fin r => (l' : ℕ) < (l : ℕ) + 1))
      = insert l (Finset.univ.filter (fun l' : Fin r => (l' : ℕ) < (l : ℕ))) := by
    ext l'
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_insert]
    constructor
    · intro h
      rcases lt_or_ge (l' : ℕ) (l : ℕ) with hlt | hge
      · exact Or.inr hlt
      · exact Or.inl (Fin.val_injective (by omega))
    · rintro (rfl | h) <;> omega
  have hnot : l ∉ (Finset.univ.filter (fun l' : Fin r => (l' : ℕ) < (l : ℕ))) := by simp
  rw [hins, Finset.sum_insert hnot]
  ring

end SingleWeight

/-! ### 2. The eigenvalue limit at the outlier index -/

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-- **The eigenvalue limit.** Under `EigSep` the sorted eigenvalue of `X_W X_Wᵀ` at the index
`l` tends in probability to the outlier `swRho c w (γ l)`. The aligned mirror is
`tendstoInProb_eigVal_het` (`RankR/Het/Align.lean:322`); the sorted index there is
`Scalars.ellSup`, here it is `l` itself, because `EigSep.sorted` lists the roots in strictly
decreasing order and `swRho` is strictly increasing in the root.

Route: at `δ' = min δ (ε / 2)` the two count events of unit G3
(`tendsto_measure_count_Ioi_tau_sw`) at `ρ_l ∓ δ'`, with `s = l + 1` and `s = l`, trap the
eigenvalue in `(ρ_l - δ', ρ_l + δ']`, so `|λ - ρ_l| ≤ δ' < ε` outside the union of their
complements. The index is in range for large `N` (`eventually_r_le_min`), so `eigVal` is the
sorted eigenvalue there. -/
theorem tendstoInProb_eigVal_sw [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)
    (hedge : m.HeteroEdgeR w c (MPhet.bHet c w))
    {γ : Fin r → ℝ} {z : Fin r → EuclideanSpace ℝ (Fin r)}
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z) (l : Fin r) :
    TendstoInProb μ
      (fun N ω => eigVal (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
        (isHermitian_mul_transpose_self (m.stackXW w N ω)) (l : ℕ))
      (SingleWeight.swRho c w (γ l)) := by
  have hρ : ∀ l' : Fin r, MPhet.bHet c w < SingleWeight.swRho c w (γ l') := fun l' =>
    SingleWeight.EigSep.bHet_lt_swRho hc hw hsep l'
  have hanti : StrictAnti fun l' => SingleWeight.swRho c w (γ l') :=
    SingleWeight.EigSep.strictAnti_swRho hc hw hsep
  obtain ⟨δ, hδ, hgapb, hgap2⟩ := SingleWeight.exists_margin_sw hρ hanti.injective
  have hrn : ∀ᶠ N in atTop, r ≤ ∑ i, n i N :=
    (m.eventually_r_le_min hreg).mono fun N h => (le_min_iff.mp h).1
  intro ε hε
  obtain ⟨δ', hδ'def⟩ : ∃ t : ℝ, t = min δ (ε / 2) := ⟨_, rfl⟩
  have hδ'pos : 0 < δ' := by rw [hδ'def]; exact lt_min hδ (by linarith)
  have hδ'le : δ' ≤ δ := by rw [hδ'def]; exact min_le_left _ _
  have hδ'ε : δ' ≤ ε / 2 := by rw [hδ'def]; exact min_le_right _ _
  have hgap2' : ∀ k l' : Fin r, k ≠ l' →
      2 * δ' < |SingleWeight.swRho c w (γ k) - SingleWeight.swRho c w (γ l')| := by
    intro k l' hkl
    have := hgap2 k l' hkl
    linarith
  obtain ⟨τm, hτmdef⟩ : ∃ t : ℝ, t = SingleWeight.swRho c w (γ l) - δ' := ⟨_, rfl⟩
  obtain ⟨τp, hτpdef⟩ : ∃ t : ℝ, t = SingleWeight.swRho c w (γ l) + δ' := ⟨_, rfl⟩
  have hτmb : MPhet.bHet c w < τm := by
    rw [hτmdef]; have := hgapb l; linarith
  have hτpb : MPhet.bHet c w < τp := by
    rw [hτpdef]; have := hgapb l; linarith
  have hgapm : ∀ l' : Fin r,
      ((l' : ℕ) < (l : ℕ) + 1 → τm + δ' / 2 ≤ SingleWeight.swRho c w (γ l'))
        ∧ ((l : ℕ) + 1 ≤ (l' : ℕ) → SingleWeight.swRho c w (γ l') + δ' / 2 < τm) := by
    rw [hτmdef]
    exact SingleWeight.gap_sub_of_margin hanti hδ'pos hgap2' l
  have hgapp : ∀ l' : Fin r,
      ((l' : ℕ) < (l : ℕ) → τp + δ' / 2 ≤ SingleWeight.swRho c w (γ l'))
        ∧ ((l : ℕ) ≤ (l' : ℕ) → SingleWeight.swRho c w (γ l') + δ' / 2 < τp) := by
    rw [hτpdef]
    exact SingleWeight.gap_add_of_margin hanti hδ'pos hgap2' l
  -- the two count events of unit G3, and their complements
  have hAm := m.tendsto_measure_count_Ioi_tau_sw w c hc hw hreg hG hpd hedge hsep
    (mg := δ' / 2) (by linarith) hτmb (s := (l : ℕ) + 1) l.isLt hgapm
  have hAp := m.tendsto_measure_count_Ioi_tau_sw w c hc hw hreg hG hpd hedge hsep
    (mg := δ' / 2) (by linarith) hτpb (s := (l : ℕ)) l.isLt.le hgapp
  have hBm : Tendsto (fun N => μ N ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
      (τm < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
        ↔ (q : ℕ) < (l : ℕ) + 1)})ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_count_Ioi_het w N τm ((l : ℕ) + 1)).nullMeasurableSet) hAm
  have hBp : Tendsto (fun N => μ N ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
      (τp < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
        ↔ (q : ℕ) < (l : ℕ))})ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_count_Ioi_het w N τp (l : ℕ)).nullMeasurableSet) hAp
  -- the union bound, on the eventual range of the index
  refine OutliersR.tendsto_measure_zero_of_eventually_subset
    (t := fun N => ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
        (τm < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
          ↔ (q : ℕ) < (l : ℕ) + 1)})ᶜ
      ∪ ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
        (τp < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
          ↔ (q : ℕ) < (l : ℕ))})ᶜ)
    ?_ (tendsto_measure_zero_union hBm hBp)
  filter_upwards [hrn] with N hN
  intro ω hω
  by_cases hlow : ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
      (τm < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
        ↔ (q : ℕ) < (l : ℕ) + 1)
  · by_cases hhigh : ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
        (τp < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
          ↔ (q : ℕ) < (l : ℕ))
    · exfalso
      have hℓ' : (l : ℕ) < Fintype.card (Fin (∑ i, n i N)) := by
        rw [Fintype.card_fin]; exact lt_of_lt_of_le l.isLt hN
      have hω' : ε ≤ |eigVal (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
          (isHermitian_mul_transpose_self (m.stackXW w N ω)) (l : ℕ)
          - SingleWeight.swRho c w (γ l)| := hω
      rw [eigVal_eq _ _ hℓ'] at hω'
      have h1 : τm < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀
          ⟨(l : ℕ), hℓ'⟩ := (hlow ⟨(l : ℕ), hℓ'⟩).mpr (Nat.lt_succ_self _)
      have h2 : ¬ τp < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀
          ⟨(l : ℕ), hℓ'⟩ := fun h => lt_irrefl _ ((hhigh ⟨(l : ℕ), hℓ'⟩).mp h)
      rw [not_lt] at h2
      rw [hτmdef] at h1
      rw [hτpdef] at h2
      have habs : |(isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀
          ⟨(l : ℕ), hℓ'⟩ - SingleWeight.swRho c w (γ l)| ≤ δ' :=
        abs_le.mpr ⟨by linarith, by linarith⟩
      linarith
    · exact Set.mem_union_right _ hhigh
  · exact Set.mem_union_left _ hlow

/-! ### 3. The target: the `align` field of `SingleWeightLaw` -/

/-- **Unit G4, the target** (the `align` field of `SingleWeightLaw`,
`RankR/SingleWeight/Main.lean:67`; the paper display is `main_paper.tex:2153` inside
`prop:gen_rank_stacksvd_singleweight`, `:2112`). Under `EigSep` the overlap of the sorted
right singular direction `l` of the weighted stack `X_W` with the spike direction `v_k` tends
in probability to `swTerm … (γ l) (z l) * (z_l)_k²`.

Route (the sandwich): the two thresholds `ρ_l ∓ δ` of `SingleWeight.exists_margin_sw`, the two
count events of unit G3 at both (`s = l + 1` and `s = l`), the eigenvalue limit
`tendstoInProb_eigVal_sw`, the two half-line limits of unit G3
(`tendstoInProb_normSq_specProj_Ioi_tau_sw`) whose difference is the single term
`(z_l)_k² / ν_l` (`SingleWeight.sum_filter_lt_succ_sub`), the quotient, and on the
intersection of the count events the duality of task E2
(`Het.overlapIdx_eq_normSq_specProjIdx_div`) with the window lemmas
`specProjIdx_eq_specProj_Ioc` and `normSq_specProj_Ioc`. The overlap identity
`1/(ρ_l ν_l) = swTerm` is `SingleWeight.one_div_swRho_mul_swNu` (unit G0). The aligned mirror
is `align_cross_het_of_gaussian_sup` (`RankR/Het/Align.lean:444`); there the limit is `0` off
the diagonal, here `EigSep` makes every component supercritical and each pair `(l, k)` carries
its own value. -/
theorem align_sw_of_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r)
    (hedge : m.HeteroEdgeR w c (MPhet.bHet c w))
    {γ : Fin r → ℝ} {z : Fin r → EuclideanSpace ℝ (Fin r)}
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z) (l k : Fin r) :
    TendstoInProb μ
      (fun N ω => overlapIdx (m.stackXW w N ω) (l : ℕ) (m.colVecG N k))
      (SingleWeight.swTerm (fun i => (m.tbl i).θ) m.R w c (γ l) (z l)
        * (WithLp.ofLp (z l) k) ^ 2) := by
  have hρ : ∀ l' : Fin r, MPhet.bHet c w < SingleWeight.swRho c w (γ l') := fun l' =>
    SingleWeight.EigSep.bHet_lt_swRho hc hw hsep l'
  have hanti : StrictAnti fun l' => SingleWeight.swRho c w (γ l') :=
    SingleWeight.EigSep.strictAnti_swRho hc hw hsep
  obtain ⟨δ, hδ, hgapb, hgap2⟩ := SingleWeight.exists_margin_sw hρ hanti.injective
  have hrn : ∀ᶠ N in atTop, r ≤ ∑ i, n i N :=
    (m.eventually_r_le_min hreg).mono fun N h => (le_min_iff.mp h).1
  obtain ⟨τm, hτmdef⟩ : ∃ t : ℝ, t = SingleWeight.swRho c w (γ l) - δ := ⟨_, rfl⟩
  obtain ⟨τp, hτpdef⟩ : ∃ t : ℝ, t = SingleWeight.swRho c w (γ l) + δ := ⟨_, rfl⟩
  have hτmb : MPhet.bHet c w < τm := by
    rw [hτmdef]; have := hgapb l; linarith
  have hτpb : MPhet.bHet c w < τp := by
    rw [hτpdef]; have := hgapb l; linarith
  have hab : τm ≤ τp := by rw [hτmdef, hτpdef]; linarith
  have hgapm : ∀ l' : Fin r,
      ((l' : ℕ) < (l : ℕ) + 1 → τm + δ / 2 ≤ SingleWeight.swRho c w (γ l'))
        ∧ ((l : ℕ) + 1 ≤ (l' : ℕ) → SingleWeight.swRho c w (γ l') + δ / 2 < τm) := by
    rw [hτmdef]
    exact SingleWeight.gap_sub_of_margin hanti hδ hgap2 l
  have hgapp : ∀ l' : Fin r,
      ((l' : ℕ) < (l : ℕ) → τp + δ / 2 ≤ SingleWeight.swRho c w (γ l'))
        ∧ ((l : ℕ) ≤ (l' : ℕ) → SingleWeight.swRho c w (γ l') + δ / 2 < τp) := by
    rw [hτpdef]
    exact SingleWeight.gap_add_of_margin hanti hδ hgap2 l
  -- step 1: the two count events of unit G3, and their complements
  have hAm := m.tendsto_measure_count_Ioi_tau_sw w c hc hw hreg hG hpd hedge hsep
    (mg := δ / 2) (by linarith) hτmb (s := (l : ℕ) + 1) l.isLt hgapm
  have hAp := m.tendsto_measure_count_Ioi_tau_sw w c hc hw hreg hG hpd hedge hsep
    (mg := δ / 2) (by linarith) hτpb (s := (l : ℕ)) l.isLt.le hgapp
  have hBm : Tendsto (fun N => μ N ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
      (τm < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
        ↔ (q : ℕ) < (l : ℕ) + 1)})ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_count_Ioi_het w N τm ((l : ℕ) + 1)).nullMeasurableSet) hAm
  have hBp : Tendsto (fun N => μ N ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
      (τp < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
        ↔ (q : ℕ) < (l : ℕ))})ᶜ) atTop (𝓝 0) :=
    tendsto_measure_compl_zero
      (fun N => (m.measurableSet_count_Ioi_het w N τp (l : ℕ)).nullMeasurableSet) hAp
  -- step 2: the two half-line limits of unit G3 and the eigenvalue limit
  have hg1 := m.tendstoInProb_normSq_specProj_Ioi_tau_sw w c hc hw hreg hG hpd hedge hsep
    (mg := δ / 2) (by linarith) hτmb (s := (l : ℕ) + 1) l.isLt hgapm k
  have hg2 := m.tendstoInProb_normSq_specProj_Ioi_tau_sw w c hc hw hreg hG hpd hedge hsep
    (mg := δ / 2) (by linarith) hτpb (s := (l : ℕ)) l.isLt.le hgapp k
  have hlam := m.tendstoInProb_eigVal_sw w c hc hw hreg hG hpd hedge hsep l
  have hρpos : 0 < SingleWeight.swRho c w (γ l) :=
    SingleWeight.swRho_pos hc hw (hsep.root l).1 (hsep.thresh l)
  have hνpos : 0 < SingleWeight.swNu (fun i => (m.tbl i).θ) m.R c w (γ l) (z l) :=
    SingleWeight.EigSep.swNu_pos hc hw hsep l
  -- step 3: the quotient limit, one term after the difference of the two half-line sums
  have hq : TendstoInProb μ
      (fun N ω => (‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ) (Set.Ioi τm)
          (WithLp.toLp 2 fun q => m.QmatHetR w N ω q k)‖ ^ 2
        - ‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ) (Set.Ioi τp)
          (WithLp.toLp 2 fun q => m.QmatHetR w N ω q k)‖ ^ 2)
        / eigVal (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
          (isHermitian_mul_transpose_self (m.stackXW w N ω)) (l : ℕ))
      (SingleWeight.swTerm (fun i => (m.tbl i).θ) m.R w c (γ l) (z l)
        * (WithLp.ofLp (z l) k) ^ 2) := by
    refine FormsR.tendstoInProb_congr_limit ?_ ((hg1.sub hg2).div hlam hρpos.ne')
    rw [SingleWeight.sum_filter_lt_succ_sub,
      ← SingleWeight.one_div_swRho_mul_swNu hc hw (hsep.root l).1 (hsep.thresh l)
        (hsep.eigvec l).1 (hsep.eigvec l).2]
    field_simp
  -- step 4: the union bound, with the duality of task E2 on the intersection
  refine TendstoInProb.of_tendsto_measure_ne_of_tendsto ?_ hq
  refine OutliersR.tendsto_measure_zero_of_eventually_subset
    (t := fun N => ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
        (τm < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
          ↔ (q : ℕ) < (l : ℕ) + 1)})ᶜ
      ∪ ({ω : Ω N | ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
        (τp < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
          ↔ (q : ℕ) < (l : ℕ))})ᶜ)
    ?_ (tendsto_measure_zero_union hBm hBp)
  filter_upwards [hrn] with N hN
  intro ω hω
  by_cases hlow : ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
      (τm < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
        ↔ (q : ℕ) < (l : ℕ) + 1)
  · by_cases hhigh : ∀ q : Fin (Fintype.card (Fin (∑ i, n i N))),
        (τp < (isHermitian_mul_transpose_self (m.stackXW w N ω)).eigenvalues₀ q
          ↔ (q : ℕ) < (l : ℕ))
    · refine absurd ?_ hω
      have hℓ' : (l : ℕ) < Fintype.card (Fin (∑ i, n i N)) := by
        rw [Fintype.card_fin]; exact lt_of_lt_of_le l.isLt hN
      have hℓ : (l : ℕ) < min (∑ i, n i N) (d N) :=
        lt_min (lt_of_lt_of_le l.isLt hN) (by have h1 := hpd N; have h2 := l.isLt; omega)
      have hpos' : 0 < eigVal (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
          (isHermitian_mul_transpose_self (m.stackXW w N ω)) (l : ℕ) := by
        rw [eigVal_eq _ _ hℓ']
        exact lt_trans (lt_trans (MPhet.bHet_pos hc hw) hτmb)
          ((hlow ⟨(l : ℕ), hℓ'⟩).mpr (Nat.lt_succ_self _))
      have hpos : 0 < eigVal ((m.stackXW w N ω)ᵀ * m.stackXW w N ω)
          (isHermitian_transpose_mul_self (m.stackXW w N ω)) (l : ℕ) := by
        rw [Het.eigVal_gram_comm (m.stackXW w N ω) hℓ]
        exact hpos'
      change overlapIdx (m.stackXW w N ω) (l : ℕ) (m.colVecG N k)
        = (‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ) (Set.Ioi τm)
            (WithLp.toLp 2 fun q => m.QmatHetR w N ω q k)‖ ^ 2
          - ‖specProj (m.stackXW w N ω * (m.stackXW w N ω)ᵀ) (Set.Ioi τp)
            (WithLp.toLp 2 fun q => m.QmatHetR w N ω q k)‖ ^ 2)
          / eigVal (m.stackXW w N ω * (m.stackXW w N ω)ᵀ)
            (isHermitian_mul_transpose_self (m.stackXW w N ω)) (l : ℕ)
      rw [Het.overlapIdx_eq_normSq_specProjIdx_div (m.stackXW w N ω) hpos (m.colVecG N k),
        m.stackXW_mulVec_colVecG w N ω k,
        specProjIdx_eq_specProj_Ioc (isHermitian_mul_transpose_self (m.stackXW w N ω)) hlow
          hhigh,
        normSq_specProj_Ioc (isHermitian_mul_transpose_self (m.stackXW w N ω)) hab]
    · exact Set.mem_union_right _ hhigh
  · exact Set.mem_union_left _ hlow

end UnalignedModelR

end StackedSVD

-- Deciding check of unit G4 (`notes/archive/trackG_plan.md`), both elaborated on 2026-09-05:
-- #check @StackedSVD.UnalignedModelR.tendstoInProb_eigVal_sw
-- #check @StackedSVD.UnalignedModelR.align_sw_of_gaussian
