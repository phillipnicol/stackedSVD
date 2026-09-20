/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.FormsBridge
import StackedSVD.RMT.General.Iso
import StackedSVD.RMT.General.ProbC
import StackedSVD.RMT.General.Defs

/-!
# The six fields of `ResolventFormsC` from six scalar limits (Stage 1, unit F1b)

`notes/archive/prop_single_table_general.md` section 5, unit F1b. The file is pure algebra. It takes
the six scalar limits of the isotropic forms on the Wishart block `gram Z` and on its
companion `gramC Z`, rewrites each of the six forms on the downdated matrix `W₀` through the
deterministic identities A1 to A7 of `RMT/General/FormsBridge.lean`, and passes to the limit
with the closure lemmas of `RMT/General/ProbC.lean`.

The scalar dictionary, at a fixed `z` with `0 < z.im`, with `W = gram Z`, `Wc = gramC Z`:

| form | limit | hypothesis |
|---|---|---|
| `qformC W z v` | `mC c z` | `hL1` |
| `qformC Wc z u` | `mTildeC c z` | `hL2` |
| `cformC W z v g` | `0` | `hL3` |
| `qform2C W z v` | `mCDeriv c z` | `hL1'` |
| `qform2C Wc z u` | `mTildeCDeriv c z` | `hL2'` |
| `cform2C W z v g` | `0` | `hL3'` |
-/

open MeasureTheory Filter Topology
open scoped Matrix

namespace StackedSVD

/-- The second complex resolvent form is symmetric. Copy of `StackedSVD.FormsR.cform2C_comm`
(`RankR/RMT/Forms.lean:236`), which lives in a rank-`r` module the general single-table route
may not import; `GenRMT.cformC_comm` (`RMT/General/Iso.lean:513`) is the first-order twin and
is reachable. -/
private theorem cform2C_comm_copy {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ} {z : ℂ}
    (hW : W.IsHermitian) (hz : z.im ≠ 0) (x y : Fin D → ℝ) :
    R4C.cform2C W z x y = R4C.cform2C W z y x := by
  rw [R4C.cform2C_eq_sum hW hz, R4C.cform2C_eq_sum hW hz]
  exact Finset.sum_congr rfl fun a _ => by ring_nf

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- **Unit F1b.** The six fields of `SpikedModel.ResolventFormsC`
(`RMT/General/Defs.lean`) from the six scalar limits on the Wishart block
`gram Z` and its companion `gramC Z`. Deterministic input: the identities A1 to A7 of
`RMT/General/FormsBridge.lean` and the scalar identities B3 to B5 of the same file. -/
theorem resolventFormsC_of_limits {c : ℝ}
    (hc : 0 < c) (m : SpikedModel μ n d)
    (hL1 : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω => ‖R4C.qformC (GenRMT.gram (m.Z N ω)) z
      (WithLp.ofLp (m.v N)) - MP.mC c z‖) 0)
    (hL2 : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω => ‖R4C.qformC (GenRMT.gramC (m.Z N ω)) z
      (WithLp.ofLp (m.u N)) - MP.mTildeC c z‖) 0)
    (hL3 : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω => ‖R4C.cformC (GenRMT.gram (m.Z N ω)) z
      (WithLp.ofLp (m.v N)) (m.gvec N ω)‖) 0)
    (hL1' : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω => ‖R4C.qform2C (GenRMT.gram (m.Z N ω)) z
      (WithLp.ofLp (m.v N)) - MP.mCDeriv c z‖) 0)
    (hL2' : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω => ‖R4C.qform2C (GenRMT.gramC (m.Z N ω)) z
      (WithLp.ofLp (m.u N)) - MP.mTildeCDeriv c z‖) 0)
    (hL3' : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω => ‖R4C.cform2C (GenRMT.gram (m.Z N ω)) z
      (WithLp.ofLp (m.v N)) (m.gvec N ω)‖) 0) :
    m.ResolventFormsC c := by
  -- the companion denominator `z b` and its nonzero limit `z m̃`
  have hzm : ∀ z : ℂ, 0 < z.im → z * MP.mTildeC c z ≠ 0 :=
    fun _ hz => MP.z_mul_mTildeC_ne_zero hc hz
  have hden : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω =>
      ‖z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))
        - z * MP.mTildeC c z‖) 0 :=
    fun z hz => tendstoInProbC_const_mul z (hL2 z hz)
  -- the two cross forms in the `- 0` shape the closure lemmas take, and their transposes
  have hC : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω =>
      ‖R4C.cformC (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N)) (m.gvec N ω) - 0‖) 0 := by
    intro z hz
    simpa only [sub_zero] using hL3 z hz
  have hCs : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω =>
      ‖R4C.cformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) (WithLp.ofLp (m.v N)) - 0‖) 0 := by
    intro z hz
    exact tendstoInProbC_congr
      (fun N ω => GenRMT.cformC_comm (GenRMT.gram_isHermitian (m.Z N ω)) hz.ne' _ _) (hC z hz)
  have hD : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω =>
      ‖R4C.cform2C (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N)) (m.gvec N ω) - 0‖) 0 := by
    intro z hz
    simpa only [sub_zero] using hL3' z hz
  have hDs : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω =>
      ‖R4C.cform2C (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) (WithLp.ofLp (m.v N)) - 0‖) 0 := by
    intro z hz
    exact tendstoInProbC_congr
      (fun N ω => cform2C_comm_copy (GenRMT.gram_isHermitian (m.Z N ω)) hz.ne' _ _) (hD z hz)
  -- the square of the companion denominator, `(z b)²`, and its limit `(z m̃)²`
  have hsq2 : ∀ z : ℂ, 0 < z.im → TendstoInProb μ (fun N ω =>
      ‖(z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))) ^ 2
        - (z * MP.mTildeC c z) ^ 2‖) 0 := by
    intro z hz
    have h := tendstoInProbC_mul (hden z hz) (hden z hz)
    rw [show z * MP.mTildeC c z * (z * MP.mTildeC c z) = (z * MP.mTildeC c z) ^ 2 from by ring]
      at h
    exact tendstoInProbC_congr (fun N ω => (pow_two _).symm) h
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  · -- **vvC**, item A1 at `x = v`: the correction is a product of two null forms over `z b`.
    intro z hz
    have h1 := tendstoInProbC_mul (hC z hz) (hCs z hz)
    have h2 := tendstoInProbC_div (hzm z hz) h1 (hden z hz)
    have h3 := tendstoInProbC_sub (hL1 z hz) h2
    have hlim : MP.mC c z - 0 * 0 / (z * MP.mTildeC c z) = MP.mC c z := by ring
    rw [hlim] at h3
    exact tendstoInProbC_congr
      (fun N ω => (m.qformC_W0_eq N ω hz (WithLp.ofLp (m.v N))).symm) h3
  · -- **ggC**, item A2 and its scalar twin B4.
    intro z hz
    have h1 := tendstoInProbC_inv (hzm z hz) (hden z hz)
    have h2 := tendstoInProbC_sub (tendstoInProbC_const μ (-1 : ℂ)) h1
    have hlim : (-1 : ℂ) - (z * MP.mTildeC c z)⁻¹ = MP.mC c z :=
      (MP.mC_eq_neg_one_sub_inv hc hz).symm
    rw [hlim] at h2
    exact tendstoInProbC_congr (fun N ω => (m.qformC_W0_gvec_eq N ω hz).symm) h2
  · -- **vgC**, item A3: a null form over `z b`.
    intro z hz
    have h1 := tendstoInProbC_div (hzm z hz) (hC z hz) (hden z hz)
    have h2 := tendstoInProbC_neg h1
    have hlim : -(0 / (z * MP.mTildeC c z)) = (0 : ℂ) := by ring
    rw [hlim] at h2
    have h3 := tendstoInProbC_congr
      (fun N ω => (m.cformC_W0_gvec_eq N ω hz (WithLp.ofLp (m.v N))).symm) h2
    simpa only [sub_zero] using h3
  · -- **vv2C**, item A5 at `x = y = v`, with `qform2C W z y = cform2C W z y y` by `rfl`.
    intro z hz
    have hA5v : ∀ (N : ℕ) (ω : Ω N),
        R4C.qform2C (m.W0 N ω) z (WithLp.ofLp (m.v N))
          = R4C.qform2C (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N))
            - (R4C.cform2C (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N)) (m.gvec N ω)
                  * R4C.cformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) (WithLp.ofLp (m.v N))
                + R4C.cformC (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N)) (m.gvec N ω)
                  * R4C.cform2C (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) (WithLp.ofLp (m.v N)))
                / (z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)))
            + R4C.cformC (GenRMT.gram (m.Z N ω)) z (WithLp.ofLp (m.v N)) (m.gvec N ω)
                * (R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))
                    + z * R4C.qform2C (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N)))
                * R4C.cformC (GenRMT.gram (m.Z N ω)) z (m.gvec N ω) (WithLp.ofLp (m.v N))
                / (z * R4C.qformC (GenRMT.gramC (m.Z N ω)) z (WithLp.ofLp (m.u N))) ^ 2 :=
      fun N ω => m.cform2C_W0_eq N ω hz (WithLp.ofLp (m.v N)) (WithLp.ofLp (m.v N))
    have n3 := tendstoInProbC_add (tendstoInProbC_mul (hD z hz) (hCs z hz))
      (tendstoInProbC_mul (hC z hz) (hDs z hz))
    have n4 := tendstoInProbC_div (hzm z hz) n3 (hden z hz)
    have s1 := tendstoInProbC_sub (hL1' z hz) n4
    have p1 := tendstoInProbC_add (hL2 z hz) (tendstoInProbC_const_mul z (hL2' z hz))
    have p3 := tendstoInProbC_mul (tendstoInProbC_mul (hC z hz) p1) (hCs z hz)
    have p5 := tendstoInProbC_div (pow_ne_zero 2 (hzm z hz)) p3 (hsq2 z hz)
    have s2 := tendstoInProbC_add s1 p5
    have hlim : MP.mCDeriv c z - (0 * 0 + 0 * 0) / (z * MP.mTildeC c z)
        + 0 * (MP.mTildeC c z + z * MP.mTildeCDeriv c z) * 0 / (z * MP.mTildeC c z) ^ 2
        = MP.mCDeriv c z := by ring
    rw [hlim] at s2
    exact tendstoInProbC_congr (fun N ω => (hA5v N ω).symm) s2
  · -- **gg2C**, item A6 and its scalar twin B5.
    intro z hz
    have p1 := tendstoInProbC_add (hL2 z hz) (tendstoInProbC_const_mul z (hL2' z hz))
    have q := tendstoInProbC_div (pow_ne_zero 2 (hzm z hz)) p1 (hsq2 z hz)
    rw [← MP.mCDeriv_eq_of_mTildeC hc hz] at q
    exact tendstoInProbC_congr (fun N ω => (m.qform2C_W0_gvec_eq N ω hz).symm) q
  · -- **vg2C**, item A7 at `x = v`: every term carries a null form.
    intro z hz
    have w1 := tendstoInProbC_add (tendstoInProbC_const μ (1 : ℂ)) (hden z hz)
    have p1 := tendstoInProbC_add (hL2 z hz) (tendstoInProbC_const_mul z (hL2' z hz))
    have t3 := tendstoInProbC_mul (hC z hz) p1
    have t4 := tendstoInProbC_add (tendstoInProbC_mul (hD z hz) w1) t3
    have t5 := tendstoInProbC_div (hzm z hz) t4 (hden z hz)
    have t6 := tendstoInProbC_sub (hD z hz) t5
    have t8 := tendstoInProbC_div (pow_ne_zero 2 (hzm z hz)) (tendstoInProbC_mul t3 w1) (hsq2 z hz)
    have t9 := tendstoInProbC_add t6 t8
    have hlim : 0 - (0 * (1 + z * MP.mTildeC c z)
          + 0 * (MP.mTildeC c z + z * MP.mTildeCDeriv c z)) / (z * MP.mTildeC c z)
        + 0 * (MP.mTildeC c z + z * MP.mTildeCDeriv c z) * (1 + z * MP.mTildeC c z)
          / (z * MP.mTildeC c z) ^ 2 = (0 : ℂ) := by ring
    rw [hlim] at t9
    have h := tendstoInProbC_congr
      (fun N ω => (m.cform2C_W0_gvec_eq N ω hz (WithLp.ofLp (m.v N))).symm) t9
    simpa only [sub_zero] using h

end SpikedModel

end StackedSVD
