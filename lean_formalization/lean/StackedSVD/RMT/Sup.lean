/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.R1
import StackedSVD.RMT.R2
import StackedSVD.RMT.R3
import StackedSVD.RMT.R5
import StackedSVD.RMT.Symmetry
import StackedSVD.RMT.Simplicity

/-!
# Item Sup: the Layer 2 assembly (milestone L2a)

This file plugs items R0, R1, R2, R3, T, Sym and S into the analytic core R5. It has two
theorems.

* `SpikedModel.resolventLimits_of_gaussian` builds the 7-field interface
  `SpikedModel.ResolventLimits` (decision D6) for a Gaussian spiked model in the proportional
  regime. The `edge` field is item R3. The six form fields come from item T at every real
  `x > bulkEdge c`, whose complex inputs at `x + i η` are the six conclusions of item R2, whose
  own trace inputs are R1a and R1b.
* `SpikedModel.singleTableLaw_of_gaussian_supercritical` is **milestone L2a** of `PLAN.md`:
  the full `SingleTableLaw` for `c < θ⁴`. It differs from `singleTableLaw_of_gaussian`
  (stated in `RMT.lean`, proved in `RMT/Full.lean` on 2026-08-30) by the two hypotheses
  `hθ : c < m.θ ^ 4` and `hn2 : ∀ N, 2 ≤ n N`.

The subcritical half is not touched here: R3⁻ and R6' close it later.

Three bridges connect the three spellings of the Wishart block.

1. `R1.W0 = R2.W0` (`R1_W0_eq_R2_W0`, `rfl`): the two files were written in parallel.
2. `R2.W0 (B N ω) = m.W0 N ω`, from `SpikedModel.exists_block_hasLaw` (item R0). Not `rfl`:
   the block `B` comes from a rotation, so the identity is propositional.
3. `R2.gOf ((Z N ω)ᵀ u) = m.gvec N ω` (`gOf_eq_gvec`), from `SpikedModel.gvec_eq_smul`.

STATUS: see `notes/archive/agent_reports/proof_sup.md`. The `sorryAx` this file inherits comes only
from the three tracked `sorry` of `RMT/R1.lean`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace ENNReal

namespace StackedSVD

/-! ### Bridge 1: the two spellings of the Wishart block agree -/

/-- `R1.W0` and `R2.W0` are the same definition (`R2` does not import `R1`). -/
theorem R1_W0_eq_R2_W0 {p D : ℕ} (Y : Matrix (Fin p) (Fin D) ℝ) : R1.W0 Y = R2.W0 Y := rfl

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-! ### Bridge 3: `R2.gOf` of the unscaled vector is the model's `g` -/

/-- `g = d^{-1/2} Zᵀ u` in the spelling item R2 uses. -/
theorem gOf_eq_gvec (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) :
    R2.gOf ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)) = m.gvec N ω :=
  (m.gvec_eq_smul N ω).symm

/-- `v ⬝ᵥ v = 1`, the model's own unit vector in dot-product form. -/
theorem dotProduct_v_self (m : SpikedModel μ n d) (N : ℕ) :
    WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) = 1 :=
  R2.dotProduct_ofLp_self (m.hv N)

/-! ### The interface, for Gaussian noise -/

/-- **The 7-field `ResolventLimits` interface for a Gaussian spiked model.** `hn2` is decision
D11: item R0 splits `n N` rows into `1 + (n N - 1)`, and R1, R2, R3 need the block to be
nonempty. -/
theorem resolventLimits_of_gaussian [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) (hG : m.GaussianNoise) (hn2 : ∀ N, 2 ≤ n N) :
    m.ResolventLimits c := by
  classical
  have hd : Tendsto d atTop atTop := hreg.2.1
  have hvdot : ∀ N, WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) = 1 := m.dotProduct_v_self
  -- 1. the Gaussian block of item R0, with `n N = pN N + 1` and `pN N ≥ 1`.
  obtain ⟨pN, hpsucc⟩ : ∃ pN : ℕ → ℕ, ∀ N, n N = pN N + 1 :=
    ⟨fun N => n N - 1, fun N => (Nat.succ_pred_eq_of_pos (m.hn N)).symm⟩
  have hppos : ∀ N, 0 < pN N := by
    intro N
    have h1 := hn2 N
    have h2 := hpsucc N
    omega
  have hpeq : ∀ N, pN N = n N - 1 := by
    intro N
    have h2 := hpsucc N
    omega
  have hcN : Tendsto (fun N => (pN N : ℝ) / (d N : ℝ)) atTop (𝓝 c) := by
    simpa only [hpeq] using m.tendsto_pred_ratio hreg
  choose B hBW hBlaw using fun N => m.exists_block_hasLaw N hG (hpsucc N)
  have hBsnd : ∀ N, HasLaw (B N) (gaussianMatrix (pN N) (d N)) (μ N) :=
    fun N => measurePreserving_snd.fun_comp_hasLaw (hBlaw N)
  set ZZ : ∀ N, Ω N → R2.NoiseSpace (pN N) (d N) :=
    fun N ω => ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N), B N ω) with hZZdef
  have hZZlaw : ∀ N, HasLaw (ZZ N) (R2.noiseLaw (pN N) (d N)) (μ N) := hBlaw
  -- 2. the two bridges, on the block that was just chosen.
  have hW0 : ∀ (N : ℕ) (ω : Ω N), R2.W0 (ZZ N ω).2 = m.W0 N ω := fun N ω => (hBW N ω).symm
  have hgv : ∀ (N : ℕ) (ω : Ω N), R2.gOf (ZZ N ω).1 = m.gvec N ω :=
    fun N ω => m.gOf_eq_gvec N ω
  -- 3. item R3, in both the `𝓝 1` shape (the field) and the `𝓝 0` shape (item T).
  have hedge1 : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ bulkEdge c + ε}) atTop (𝓝 1) := by
    intro ε hε
    exact R3.tendsto_measure_lamMax_le μ hc B m.W0 hBW m.isHermitian_W0 hBsnd hppos m.hd hd
      hcN hε
  have hedge : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ bulkEdge c + ε}ᶜ) atTop (𝓝 0) :=
    fun ε hε => tendsto_measure_compl_zero
      (fun N => (m.measurableSet_lamMax_W0_le N _).nullMeasurableSet) (hedge1 ε hε)
  -- 4. the six complex forms of item R2, moved onto `m.W0` and `m.gvec`.
  have hforms : ∀ z : ℂ, 0 < z.im →
      TendstoInProb μ (fun N ω =>
        ‖R4C.qformC (m.W0 N ω) z (WithLp.ofLp (m.v N)) - MP.mC c z‖) 0 ∧
      TendstoInProb μ (fun N ω =>
        ‖R4C.qformC (m.W0 N ω) z (m.gvec N ω) - MP.mC c z‖) 0 ∧
      TendstoInProb μ (fun N ω =>
        ‖R4C.cformC (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω)‖) 0 ∧
      TendstoInProb μ (fun N ω =>
        ‖R4C.qform2C (m.W0 N ω) z (WithLp.ofLp (m.v N)) - MP.mCDeriv c z‖) 0 ∧
      TendstoInProb μ (fun N ω =>
        ‖R4C.qform2C (m.W0 N ω) z (m.gvec N ω) - MP.mCDeriv c z‖) 0 ∧
      TendstoInProb μ (fun N ω =>
        ‖R4C.cform2C (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω)‖) 0 := by
    intro z hz
    have hR1a : TendstoInProb μ
        (fun N ω => ‖R4C.stieltjesC (R2.W0 (ZZ N ω).2) z - MP.mC c z‖) 0 :=
      R1.tendstoInProb_stieltjesC hc hz hd hppos hcN (fun N ω => (ZZ N ω).2) hBsnd
    have hR1b : TendstoInProb μ
        (fun N ω => ‖R4C.stieltjes2C (R2.W0 (ZZ N ω).2) z - MP.mCDeriv c z‖) 0 :=
      R1.tendstoInProb_stieltjes2C hc hz hd hppos hcN (fun N ω => (ZZ N ω).2) hBsnd
    have h := R2.tendsto_forms hc hz hd hppos hcN m.v m.hv ZZ hZZlaw hR1a hR1b
    simpa only [hW0, hgv] using h
  -- 5. the norm events item T conditions on: `v ⬝ᵥ v = 1` and `g ⬝ᵥ g → 1`.
  have hgdot : TendstoInProb μ (fun N ω => m.gvec N ω ⬝ᵥ m.gvec N ω) 1 := by
    have h := R2.tendstoInProb_dotProduct_gOf ZZ hZZlaw hd
    simpa only [hgv] using h
  have hgbad :
      Tendsto (fun N => μ N {ω : Ω N | m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4}ᶜ) atTop (𝓝 0) := by
    refine tendsto_measure_zero_of_subset
      (t := fun N => {ω : Ω N | (3 : ℝ) ≤ |m.gvec N ω ⬝ᵥ m.gvec N ω - 1|})
      (fun N ω hω => ?_) (hgdot 3 (by norm_num))
    have hnot : ¬ (m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4) := hω
    have h4 : (4 : ℝ) < m.gvec N ω ⬝ᵥ m.gvec N ω := not_le.mp hnot
    change (3 : ℝ) ≤ |m.gvec N ω ⬝ᵥ m.gvec N ω - 1|
    rw [abs_of_nonneg (by linarith)]
    linarith
  have hnormvv : Tendsto (fun N => μ N {ω : Ω N |
      WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) ≤ 4 ∧
        WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) ≤ 4}ᶜ) atTop (𝓝 0) := by
    have hset : ∀ N : ℕ, {ω : Ω N |
        WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) ≤ 4 ∧
          WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) ≤ 4}ᶜ = (∅ : Set (Ω N)) := by
      intro N
      ext ω
      simp [hvdot N]
    simp only [hset, measure_empty]
    exact tendsto_const_nhds
  have hnormgg : Tendsto (fun N => μ N {ω : Ω N |
      m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4 ∧ m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4}ᶜ) atTop (𝓝 0) := by
    refine tendsto_measure_zero_of_subset (fun N ω hω => ?_) hgbad
    have hnot : ¬ (m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4 ∧ m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4) := hω
    exact fun hg => hnot ⟨hg, hg⟩
  have hnormvg : Tendsto (fun N => μ N {ω : Ω N |
      WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) ≤ 4 ∧
        m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4}ᶜ) atTop (𝓝 0) := by
    refine tendsto_measure_zero_of_subset (fun N ω hω => ?_) hgbad
    have hnot : ¬ (WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) ≤ 4 ∧
      m.gvec N ω ⬝ᵥ m.gvec N ω ≤ 4) := hω
    refine fun hg => hnot ⟨?_, hg⟩
    rw [hvdot N]
    norm_num
  have hzeroL : ∀ x : ℝ, Tendsto
      (fun η : ℝ => (fun _ : ℂ => (0 : ℂ)) ((x : ℂ) + (η : ℂ) * Complex.I)) (𝓝[>] (0 : ℝ))
      (𝓝 (((0 : ℝ) : ℂ))) := by
    intro x
    simp only [Complex.ofReal_zero]
    exact tendsto_const_nhds
  -- 6. item T, once per field.
  refine ⟨hedge1, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro x hx
    exact T.tendstoInProb_cform_of_complex' hc m.W0 m.isHermitian_W0
      (fun N _ => WithLp.ofLp (m.v N)) (fun N _ => WithLp.ofLp (m.v N)) hedge hnormvv
      (MP.mC c) (MP.m c x) (fun z hz => (hforms z hz).1) hx (MP.tendsto_mC hc hx)
  · intro x hx
    exact T.tendstoInProb_cform_of_complex' hc m.W0 m.isHermitian_W0 m.gvec m.gvec
      hedge hnormgg (MP.mC c) (MP.m c x) (fun z hz => (hforms z hz).2.1) hx
      (MP.tendsto_mC hc hx)
  · intro x hx
    exact T.tendstoInProb_cform_of_complex' hc m.W0 m.isHermitian_W0
      (fun N _ => WithLp.ofLp (m.v N)) m.gvec hedge hnormvg (fun _ => (0 : ℂ)) 0
      (fun z hz => by simpa only [sub_zero] using (hforms z hz).2.2.1) hx (hzeroL x)
  · intro x hx
    exact T.tendstoInProb_cform2_of_complex' hc m.W0 m.isHermitian_W0
      (fun N _ => WithLp.ofLp (m.v N)) (fun N _ => WithLp.ofLp (m.v N)) hedge hnormvv
      (MP.mCDeriv c) (MP.mDeriv c x) (fun z hz => (hforms z hz).2.2.2.1) hx
      (MP.tendsto_mCDeriv hc hx)
  · intro x hx
    exact T.tendstoInProb_cform2_of_complex' hc m.W0 m.isHermitian_W0 m.gvec m.gvec
      hedge hnormgg (MP.mCDeriv c) (MP.mDeriv c x) (fun z hz => (hforms z hz).2.2.2.2.1) hx
      (MP.tendsto_mCDeriv hc hx)
  · intro x hx
    exact T.tendstoInProb_cform2_of_complex' hc m.W0 m.isHermitian_W0
      (fun N _ => WithLp.ofLp (m.v N)) m.gvec hedge hnormvg (fun _ => (0 : ℂ)) 0
      (fun z hz => by simpa only [sub_zero] using (hforms z hz).2.2.2.2.2) hx (hzeroL x)

/-! ### Milestone L2a -/

/-- **Milestone L2a of `PLAN.md`.** `prop:single_table` for a Gaussian spiked model above the
threshold. Against `singleTableLaw_of_gaussian` (stated in `RMT.lean`, proved in
`RMT/Full.lean` on 2026-08-30) this adds two hypotheses: `hθ : c < m.θ ^ 4` (the supercritical
case; R3⁻ and R6' close `θ⁴ ≤ c`) and
`hn2 : ∀ N, 2 ≤ n N` (decision D11; item R0 splits one row off `n N`).

The four fields come from item R5 (`align`, `lamMax`), item Sym (`delocUniform`) and item S
(`topSimple`), on the interface `resolventLimits_of_gaussian` builds. -/
theorem singleTableLaw_of_gaussian_supercritical [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ}
    (hc : 0 < c) (m : SpikedModel μ n d) (hreg : m.Regime c) (hG : m.GaussianNoise)
    (hn2 : ∀ N, 2 ≤ n N) (hθ : c < m.θ ^ 4) :
    m.SingleTableLaw c := by
  have RL := m.resolventLimits_of_gaussian hc hreg hG hn2
  exact ⟨align_tendstoInProb RL hc hθ, delocUniform_of_gaussian m hG hreg.2.1,
    lamMax_tendstoInProb RL hc hθ, singleTableLaw_topSimple_of_gaussian m hG⟩

end SpikedModel

end StackedSVD
