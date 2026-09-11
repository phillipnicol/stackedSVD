/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.R0
import StackedSVD.RMT.R4C
import StackedSVD.RMT.MP7
import StackedSVD.Prob.NoiseLaw
import StackedSVD.Prob.NoiseMoments
import StackedSVD.Prob.LinFormMoments
import StackedSVD.Prob.PolynomialNull
import StackedSVD.Prob.Chebyshev
import StackedSVD.Prob.TendstoInProb

/-!
# Definitions and the two scalar limits, general noise law (Stage 1, unit G0)

`notes/archive/prop_single_table_general.md` section 5, unit G0. This file starts the general-law
twin of the Gaussian RMT tree under `RMT/General/`. It must not import `RMT/R1.lean`,
`RMT/R2.lean`, `RMT/SteinStep.lean`, any `Vendor/COLT83/` file, `ThetaEst.lean`, or any
other file of the Gaussian chain beyond `RMT/R0`, `RMT/R4C` and `RMT/MP7` (choice 8 of
the note); `Prob/Chebyshev.lean` supplies the two Chebyshev lemmas that `RMT/R2.lean` and
`ThetaEst.lean` each keep a private copy of.

## Content

1. `R4C.cmat'`, the rectangular twin of the square `R4C.cmat`.
2. `GenRMT.gram`, `GenRMT.gramC`, the Wishart and companion blocks at a general shape, and
   their Hermitian and rescaling facts.
3. `measurePreserving_transpose_noiseMatrix`, the general-law twin of `measurePreserving_transpose`
   (`Prob/GaussianMatrix.lean:203`): the transpose of a matrix with i.i.d. entries of law `ν`
   again has i.i.d. entries of law `ν`. The route copies `measurePreserving_transpose`'s own
   proof, through the generic `measurePreserving_uncurry` (`Prob/PolynomialNull.lean:243`) and
   a general-law twin of `measurePreserving_piCongrEquiv` (`Prob/GaussianMatrix.lean:96`).
4. Two bridges from the split objects of `RMT/R0.lean` to `GenRMT.gram`: `W0_eq_gram` writes
   `W₀` as the Gram matrix of the unscaled, downdated noise; `transpose_mul_E_eq_gram` writes
   `EᵀE` as the Gram matrix of the unscaled noise `Z`.
5. The three scalar limits that a general noise law shares with the Gaussian case:
   `tendstoInProb_dotProduct_gvec_general` (`g ⬝ g → 1`), `tendstoInProb_dotProduct_v_gvec_general`
   (`v ⬝ g → 0`) and `tendstoInProb_dotProduct_qvec_general` (`q ⬝ q → θ² + 1`). Each is the
   general-law twin of the matching Gaussian theorem of `RMT/R6.lean`. The route for the first
   two is Chebyshev on the fourth-moment bounds of `Prob/NoiseMoments.lean`, transferred from
   the canonical noise measure to the model's probability space along the `HasLaw` of
   `SpikedModel.GeneralNoise`; the third combines the first two exactly as the Gaussian file
   does, using `dotProduct_v_self` in place of the Gaussian file's rotation-invariance fact.
6. `SpikedModel.ResolventFormsC`, the six isotropic resolvent forms of `ResolventLimits` at
   complex `z`. Follow-up item F35 (2026-09-09) moved it here from
   `RMT/General/Statements.lean`. That file held the 12 target statements of the plan and
   is now retired.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace R4C

/-- A rectangular real matrix read over `ℂ`; the square case is `cmat`. -/
noncomputable def cmat' {p d : ℕ} (Y : Matrix (Fin p) (Fin d) ℝ) : Matrix (Fin p) (Fin d) ℂ :=
  Y.map (fun a => (a : ℂ))

@[simp]
theorem cmat'_apply {p d : ℕ} (Y : Matrix (Fin p) (Fin d) ℝ) (i : Fin p) (j : Fin d) :
    cmat' Y i j = (Y i j : ℂ) := rfl

/-- `cmat'` commutes with transpose. -/
theorem cmat'_transpose {p d : ℕ} (Y : Matrix (Fin p) (Fin d) ℝ) :
    cmat' Yᵀ = (cmat' Y)ᵀ := by
  ext i j
  rfl

/-- `cmat'` is `cmat` at a square shape. -/
theorem cmat'_eq_cmat {d : ℕ} (W : Matrix (Fin d) (Fin d) ℝ) : cmat' W = cmat W := rfl

end R4C

namespace GenRMT

/-- The Wishart block `d⁻¹ YᵀY`. -/
noncomputable def gram {p d : ℕ} (Y : Matrix (Fin p) (Fin d) ℝ) : Matrix (Fin d) (Fin d) ℝ :=
  ((d : ℝ))⁻¹ • (Yᵀ * Y)

/-- The companion block `d⁻¹ YYᵀ`, normalized by the same `d`. -/
noncomputable def gramC {p d : ℕ} (Y : Matrix (Fin p) (Fin d) ℝ) : Matrix (Fin p) (Fin p) ℝ :=
  ((d : ℝ))⁻¹ • (Y * Yᵀ)

theorem gram_isHermitian {p d : ℕ} (Y : Matrix (Fin p) (Fin d) ℝ) : (gram Y).IsHermitian := by
  change (((d : ℝ))⁻¹ • (Yᵀ * Y)).IsHermitian
  exact isHermitian_smul (isHermitian_transpose_mul_self Y) _

theorem gramC_isHermitian {p d : ℕ} (Y : Matrix (Fin p) (Fin d) ℝ) :
    (gramC Y).IsHermitian := by
  change (((d : ℝ))⁻¹ • (Y * Yᵀ)).IsHermitian
  exact isHermitian_smul (isHermitian_mul_transpose_self Y) _

/-- `gram Yᵀ` normalizes by `p`, not `d`, so the companion is a rescaled `gram` of the
transpose. -/
theorem gramC_eq_smul_gram_transpose {p d : ℕ} (Y : Matrix (Fin p) (Fin d) ℝ) (hp : 0 < p) :
    gramC Y = ((p : ℝ) / d) • gram Yᵀ := by
  have hpR : (p : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hp.ne'
  change (((d : ℝ))⁻¹ • (Y * Yᵀ)) = ((p : ℝ) / d) • (((p : ℝ))⁻¹ • ((Yᵀ)ᵀ * Yᵀ))
  rw [Matrix.transpose_transpose, smul_smul]
  congr 1
  rcases eq_or_ne (d : ℝ) 0 with hd0 | hd0
  · simp [hd0]
  · rw [div_eq_mul_inv, mul_comm (p : ℝ) (d : ℝ)⁻¹, mul_assoc, mul_inv_cancel₀ hpR, mul_one]

end GenRMT

/-! ### Measurability of the `mulVec`/`dotProduct` kernels

`ThetaEst.lean` proves the matching facts (`measurable_mulVec_pair`, `measurable_dotProduct_pair`,
`measurable_dotProduct_left_pair`) jointly in a varying direction, since its consumer needs
that. Here the direction is a fixed vector, so the statement is simpler; the route (unfold to
nested finite sums, `Finset.measurable_sum` and `measurable_pi_apply`) is the same. -/

private theorem measurable_mulVec {r D : ℕ} (a : Fin D → ℝ) :
    Measurable (fun Z : Matrix (Fin r) (Fin D) ℝ => Z *ᵥ a) := by
  refine measurable_pi_lambda _ fun k => ?_
  have hsum : (fun Z : Matrix (Fin r) (Fin D) ℝ => (Z *ᵥ a) k) = fun Z => ∑ l, Z k l * a l :=
    rfl
  rw [hsum]
  have hterm : ∀ l : Fin D, Measurable (fun Z : Matrix (Fin r) (Fin D) ℝ => Z k l * a l) := by
    intro l
    have h1 : Measurable (fun W : Fin r → Fin D → ℝ => W k) := measurable_pi_apply k
    have h2 : Measurable (fun v : Fin D → ℝ => v l) := measurable_pi_apply l
    have h3 : Measurable (fun W : Fin r → Fin D → ℝ => W k l) := h2.comp h1
    exact h3.mul measurable_const
  exact Finset.measurable_sum _ fun l _ => hterm l

private theorem measurable_dotProduct_mulVec {r D : ℕ} (x : Fin r → ℝ) (a : Fin D → ℝ) :
    Measurable (fun Z : Matrix (Fin r) (Fin D) ℝ => x ⬝ᵥ (Z *ᵥ a)) := by
  have hsum : (fun Z : Matrix (Fin r) (Fin D) ℝ => x ⬝ᵥ (Z *ᵥ a))
      = fun Z => ∑ k, x k * (Z *ᵥ a) k := rfl
  rw [hsum]
  have hterm : ∀ k : Fin r,
      Measurable (fun Z : Matrix (Fin r) (Fin D) ℝ => x k * (Z *ᵥ a) k) :=
    fun k => measurable_const.mul ((measurable_pi_apply k).comp (measurable_mulVec a))
  exact Finset.measurable_sum _ fun k _ => hterm k

private theorem measurable_normSqMulVec {r D : ℕ} (a : Fin D → ℝ) :
    Measurable (fun Z : Matrix (Fin r) (Fin D) ℝ => (Z *ᵥ a) ⬝ᵥ (Z *ᵥ a)) := by
  have hsum : (fun Z : Matrix (Fin r) (Fin D) ℝ => (Z *ᵥ a) ⬝ᵥ (Z *ᵥ a))
      = fun Z => ∑ k, (Z *ᵥ a) k * (Z *ᵥ a) k := rfl
  rw [hsum]
  have hterm : ∀ k : Fin r,
      Measurable (fun Z : Matrix (Fin r) (Fin D) ℝ => (Z *ᵥ a) k * (Z *ᵥ a) k) :=
    fun k => ((measurable_pi_apply k).comp (measurable_mulVec a)).mul
      ((measurable_pi_apply k).comp (measurable_mulVec a))
  exact Finset.measurable_sum _ fun k _ => hterm k

/-! ### The transpose law of `noiseMatrix` -/

/-- General-law twin of `measurePreserving_piCongrEquiv` (`Prob/GaussianMatrix.lean:96`): the
same proof, at any sigma-finite measure in place of the standard Gaussian (it only ever uses
`Measure.pi_pi`, never the specific law). -/
private theorem measurePreserving_piCongrEquiv_general {ι κ : Type*} [Fintype ι] [Fintype κ]
    (ν : Measure ℝ) [SigmaFinite ν] (e : κ ≃ ι) :
    MeasurePreserving (piCongrEquiv e) (Measure.pi fun _ : ι => ν)
      (Measure.pi fun _ : κ => ν) := by
  refine ⟨(piCongrEquiv e).measurable, ?_⟩
  refine (Measure.pi_eq fun s hs => ?_).symm
  rw [Measure.map_apply (piCongrEquiv e).measurable (MeasurableSet.univ_pi hs)]
  have hpre : (piCongrEquiv e) ⁻¹' Set.univ.pi s = Set.univ.pi fun i => s (e.symm i) := by
    ext x
    simp only [Set.mem_preimage, Set.mem_univ_pi, piCongrEquiv_apply]
    constructor
    · intro h i
      have hi := h (e.symm i)
      rwa [Equiv.apply_symm_apply] at hi
    · intro h k
      have hk := h (e k)
      rwa [Equiv.symm_apply_apply] at hk
  rw [hpre, Measure.pi_pi]
  exact Equiv.prod_comp e.symm fun k => ν (s k)

/-- General-law twin of `measurePreserving_matrixUncurry` (`Prob/GaussianMatrix.lean:60`):
`measurePreserving_uncurry` is already generic in the measure. -/
private theorem measurePreserving_matrixUncurry_general (ν : Measure ℝ) [SigmaFinite ν]
    (p d : ℕ) : MeasurePreserving (matrixUncurry p d) (noiseMatrix ν p d)
      (Measure.pi fun _ : Fin p × Fin d => ν) :=
  measurePreserving_uncurry ν

/-- **The transpose of a matrix with i.i.d. entries of a general law again has i.i.d. entries
of that law.** General-law twin of `measurePreserving_transpose` (`Prob/GaussianMatrix.lean:203`);
the proof is the same three-step composition, generalized from `gaussianReal 0 1` to `ν`. -/
theorem measurePreserving_transpose_noiseMatrix (ν : Measure ℝ) [SigmaFinite ν] (p d : ℕ) :
    MeasurePreserving (Matrix.transpose : Matrix (Fin p) (Fin d) ℝ → Matrix (Fin d) (Fin p) ℝ)
      (noiseMatrix ν p d) (noiseMatrix ν d p) := by
  have h1 := measurePreserving_matrixUncurry_general ν p d
  have h2 := measurePreserving_piCongrEquiv_general ν (Equiv.prodComm (Fin d) (Fin p))
  have h3 := (measurePreserving_matrixUncurry_general ν d p).symm (matrixUncurry d p)
  have hcomp := h3.comp (h2.comp h1)
  have hfun : ((matrixUncurry d p).symm ∘ (piCongrEquiv (Equiv.prodComm (Fin d) (Fin p))) ∘
      (matrixUncurry p d)) = (Matrix.transpose : Matrix (Fin p) (Fin d) ℝ → _) := rfl
  rwa [hfun] at hcomp

/-- Corollary in `Measure.map` form. -/
theorem noiseMatrix_map_transpose (ν : Measure ℝ) [SigmaFinite ν] (p d : ℕ) :
    (noiseMatrix ν p d).map
        (Matrix.transpose : Matrix (Fin p) (Fin d) ℝ → Matrix (Fin d) (Fin p) ℝ)
      = noiseMatrix ν d p :=
  (measurePreserving_transpose_noiseMatrix ν p d).map_eq

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- `v ⬝ᵥ v = 1`, the companion of `dotProduct_u_self` (`RMT/R0.lean:164`) for `v`. -/
theorem dotProduct_v_self (m : SpikedModel μ n d) (N : ℕ) :
    WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) = 1 := by
  have h : ⟪m.v N, m.v N⟫_ℝ = WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) :=
    inner_euclidean_eq_dotProduct (m.v N) (m.v N)
  rw [real_inner_self_eq_norm_sq, m.hv N] at h
  simpa using h.symm

/-- `W₀ = gram (√d • E⊥)`: the Wishart block is the normalized Gram matrix of the unscaled,
downdated noise. -/
theorem W0_eq_gram (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) :
    m.W0 N ω = GenRMT.gram (Real.sqrt (d N) • m.Eperp N ω) := by
  have hne : (d N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (m.hd N).ne'
  have hscal : (d N : ℝ)⁻¹ * Real.sqrt (d N) * Real.sqrt (d N) = 1 := by
    rw [mul_assoc, Real.mul_self_sqrt (Nat.cast_nonneg (d N)),
      mul_comm ((d N : ℝ))⁻¹ (d N : ℝ)]
    exact mul_inv_cancel₀ hne
  change (m.Eperp N ω)ᵀ * m.Eperp N ω
      = ((d N : ℝ))⁻¹ • ((Real.sqrt (d N) • m.Eperp N ω)ᵀ * (Real.sqrt (d N) • m.Eperp N ω))
  rw [Matrix.transpose_smul, Matrix.smul_mul, Matrix.mul_smul, smul_smul, smul_smul, hscal,
    one_smul]

/-- `EᵀE = gram Z`: the unscaled `EᵀE` is the Gram matrix of the unscaled noise. -/
theorem transpose_mul_E_eq_gram (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) :
    (m.E N ω)ᵀ * m.E N ω = GenRMT.gram (m.Z N ω) := by
  have hdd : (Real.sqrt (d N))⁻¹ * (Real.sqrt (d N))⁻¹ = ((d N : ℝ))⁻¹ := by
    rw [← mul_inv, Real.mul_self_sqrt (Nat.cast_nonneg (d N))]
  change ((Real.sqrt (d N))⁻¹ • m.Z N ω)ᵀ * ((Real.sqrt (d N))⁻¹ • m.Z N ω)
      = ((d N : ℝ))⁻¹ • ((m.Z N ω)ᵀ * m.Z N ω)
  rw [Matrix.transpose_smul, Matrix.smul_mul, Matrix.mul_smul, smul_smul, hdd]

/-! ### The three scalar limits -/

/-- `g ⬝ g → 1`. General-law twin of `tendstoInProb_dotProduct_gvec` (`RMT/R6.lean:752`). -/
theorem tendstoInProb_dotProduct_gvec_general (m : SpikedModel μ n d) {ν : Measure ℝ}
    (hν : NoiseLaw ν) (hG : m.GeneralNoise ν) (hd : Tendsto d atTop atTop) :
    TendstoInProb μ (fun N ω => m.gvec N ω ⬝ᵥ m.gvec N ω) 1 := by
  have := hν.prob
  have hK : Tendsto (fun N => ((∫ x, x ^ 4 ∂ν) + 2) / (d N : ℝ)) atTop (𝓝 0) := by
    have h1 : Tendsto (fun N => ((d N : ℝ))⁻¹) atTop (𝓝 0) :=
      (tendsto_natCast_atTop_atTop.comp hd).inv_tendsto_atTop
    have h2 := h1.const_mul ((∫ x, x ^ 4 ∂ν) + 2)
    simpa [div_eq_mul_inv, mul_comm] using h2
  have hbnd : ∀ ε : ℝ, 0 < ε → ∀ N,
      μ N {ω | ε ≤ |m.gvec N ω ⬝ᵥ m.gvec N ω - 1|}
        ≤ ENNReal.ofReal (((∫ x, x ^ 4 ∂ν) + 2) / (d N : ℝ) / ε ^ 2) := by
    intro ε hε N
    have hau : WithLp.ofLp (m.u N) ⬝ᵥ WithLp.ofLp (m.u N) = 1 := m.dotProduct_u_self N
    have hnR : (n N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (m.hn N).ne'
    have hdR : (d N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (m.hd N).ne'
    have hdd : (Real.sqrt (d N))⁻¹ * (Real.sqrt (d N))⁻¹ = ((d N : ℝ))⁻¹ := by
      rw [← mul_inv, Real.mul_self_sqrt (Nat.cast_nonneg (d N))]
    have hFmeas : Measurable (fun W : Matrix (Fin (d N)) (Fin (n N)) ℝ =>
        (d N : ℝ)⁻¹ * ((W *ᵥ WithLp.ofLp (m.u N)) ⬝ᵥ (W *ᵥ WithLp.ofLp (m.u N))) - 1) :=
      (measurable_const.mul (measurable_normSqMulVec (WithLp.ofLp (m.u N)))).sub
        measurable_const
    have hGbound := NoiseLaw.integral_sq_normSqMulVec_le (ν := ν) (r := d N) (D := n N)
      hν hau
    have hGint := NoiseLaw.integrable_sq_normSqMulVec (ν := ν) (r := d N) (D := n N) hν hau
    have hrescale : ∀ W : Matrix (Fin (d N)) (Fin (n N)) ℝ,
        ((d N : ℝ)⁻¹ * ((W *ᵥ WithLp.ofLp (m.u N)) ⬝ᵥ (W *ᵥ WithLp.ofLp (m.u N))) - 1)
          = ((n N : ℝ) / d N) *
            ((n N : ℝ)⁻¹ * ((W *ᵥ WithLp.ofLp (m.u N)) ⬝ᵥ (W *ᵥ WithLp.ofLp (m.u N)))
              - (d N : ℝ) / n N) := by
      intro W
      field_simp
    have heq : (fun W : Matrix (Fin (d N)) (Fin (n N)) ℝ =>
        ((d N : ℝ)⁻¹ * ((W *ᵥ WithLp.ofLp (m.u N)) ⬝ᵥ (W *ᵥ WithLp.ofLp (m.u N))) - 1) ^ 2)
        = fun W => ((n N : ℝ) / d N) ^ 2 *
            ((n N : ℝ)⁻¹ * ((W *ᵥ WithLp.ofLp (m.u N)) ⬝ᵥ (W *ᵥ WithLp.ofLp (m.u N)))
              - (d N : ℝ) / n N) ^ 2 := by
      funext W
      rw [hrescale W, mul_pow]
    have hFint : Integrable (fun W : Matrix (Fin (d N)) (Fin (n N)) ℝ =>
        ((d N : ℝ)⁻¹ * ((W *ᵥ WithLp.ofLp (m.u N)) ⬝ᵥ (W *ᵥ WithLp.ofLp (m.u N))) - 1) ^ 2)
        (noiseMatrix ν (d N) (n N)) := by
      rw [heq]
      exact hGint.const_mul _
    have hFbound : ∫ W, ((d N : ℝ)⁻¹
          * ((W *ᵥ WithLp.ofLp (m.u N)) ⬝ᵥ (W *ᵥ WithLp.ofLp (m.u N))) - 1) ^ 2
        ∂(noiseMatrix ν (d N) (n N)) ≤ ((∫ x, x ^ 4 ∂ν) + 2) / (d N : ℝ) := by
      rw [heq, integral_const_mul]
      have hsc : (0 : ℝ) ≤ ((n N : ℝ) / d N) ^ 2 := sq_nonneg _
      have hstep := mul_le_mul_of_nonneg_left hGbound hsc
      have hscalar : ((n N : ℝ) / d N) ^ 2 *
          (((∫ x, x ^ 4 ∂ν) + 2) * (d N) / (n N : ℝ) ^ 2)
          = ((∫ x, x ^ 4 ∂ν) + 2) / (d N : ℝ) := by
        field_simp
      rwa [hscalar] at hstep
    have hFabsmeas := hFmeas.abs
    have hFabsint : Integrable (fun W : Matrix (Fin (d N)) (Fin (n N)) ℝ =>
        |(d N : ℝ)⁻¹ * ((W *ᵥ WithLp.ofLp (m.u N)) ⬝ᵥ (W *ᵥ WithLp.ofLp (m.u N))) - 1| ^ 2)
        (noiseMatrix ν (d N) (n N)) := by
      simpa only [sq_abs] using hFint
    have hFabsC : ∫ W, |(d N : ℝ)⁻¹
          * ((W *ᵥ WithLp.ofLp (m.u N)) ⬝ᵥ (W *ᵥ WithLp.ofLp (m.u N))) - 1| ^ 2
        ∂(noiseMatrix ν (d N) (n N)) ≤ ((∫ x, x ^ 4 ∂ν) + 2) / (d N : ℝ) := by
      simpa only [sq_abs] using hFbound
    have hset : MeasurableSet {W : Matrix (Fin (d N)) (Fin (n N)) ℝ |
        ε ≤ |(d N : ℝ)⁻¹ * ((W *ᵥ WithLp.ofLp (m.u N)) ⬝ᵥ (W *ᵥ WithLp.ofLp (m.u N))) - 1|} :=
      measurableSet_le measurable_const hFabsmeas
    have hchb := Cheb.meas_ge_le_of_integral_sq _ hFabsmeas hFabsint hε hFabsC
    have hcomp := (measurePreserving_transpose_noiseMatrix ν (n N) (d N)).fun_comp_hasLaw
      (hG N)
    have hmeq := hcomp.measure_eq hset
    have hpt : ∀ ω, |(d N : ℝ)⁻¹ *
          (((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)) ⬝ᵥ ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N))) - 1|
        = |m.gvec N ω ⬝ᵥ m.gvec N ω - 1| := by
      intro ω
      have hgg : m.gvec N ω ⬝ᵥ m.gvec N ω
          = (d N : ℝ)⁻¹ *
            (((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)) ⬝ᵥ ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N))) := by
        rw [m.gvec_eq_smul N ω, smul_dotProduct, dotProduct_smul, smul_eq_mul, smul_eq_mul,
          ← mul_assoc, hdd]
      rw [hgg]
    have hseteq : {ω | ε ≤ |m.gvec N ω ⬝ᵥ m.gvec N ω - 1|}
        = {ω | ε ≤ |(d N : ℝ)⁻¹ *
            (((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)) ⬝ᵥ ((m.Z N ω)ᵀ *ᵥ WithLp.ofLp (m.u N))) - 1|} := by
      ext ω
      simp only [Set.mem_ofPred_eq, hpt]
    rw [hseteq, hmeq]
    exact hchb
  have h0 := Cheb.tendstoInProb_of_meas_le hK hbnd
  have hconv := h0.add (TendstoInProb.const μ 1)
  have hz : (0 : ℝ) + 1 = 1 := zero_add 1
  rw [hz] at hconv
  have hfun2 : ∀ (N : ℕ) (ω : Ω N),
      m.gvec N ω ⬝ᵥ m.gvec N ω - 1 + 1 = m.gvec N ω ⬝ᵥ m.gvec N ω := fun N ω => by ring
  simpa only [hfun2] using hconv

/-- `v ⬝ g → 0`. General-law twin of `tendstoInProb_dotProduct_v_gvec` (`RMT/R6.lean:767`). -/
theorem tendstoInProb_dotProduct_v_gvec_general (m : SpikedModel μ n d) {ν : Measure ℝ}
    (hν : NoiseLaw ν) (hG : m.GeneralNoise ν) (hd : Tendsto d atTop atTop) :
    TendstoInProb μ (fun N ω => WithLp.ofLp (m.v N) ⬝ᵥ m.gvec N ω) 0 := by
  have := hν.prob
  have hK : Tendsto (fun N => (1 : ℝ) / (d N : ℝ)) atTop (𝓝 0) := by
    have h1 : Tendsto (fun N => ((d N : ℝ))⁻¹) atTop (𝓝 0) :=
      (tendsto_natCast_atTop_atTop.comp hd).inv_tendsto_atTop
    simpa [one_div] using h1
  have hbnd : ∀ ε : ℝ, 0 < ε → ∀ N,
      μ N {ω | ε ≤ |WithLp.ofLp (m.v N) ⬝ᵥ m.gvec N ω - 0|}
        ≤ ENNReal.ofReal ((1 / (d N : ℝ)) / ε ^ 2) := by
    intro ε hε N
    have hau : WithLp.ofLp (m.u N) ⬝ᵥ WithLp.ofLp (m.u N) = 1 := m.dotProduct_u_self N
    have hav : WithLp.ofLp (m.v N) ⬝ᵥ WithLp.ofLp (m.v N) = 1 := m.dotProduct_v_self N
    have hdd : (Real.sqrt (d N))⁻¹ * (Real.sqrt (d N))⁻¹ = ((d N : ℝ))⁻¹ := by
      rw [← mul_inv, Real.mul_self_sqrt (Nat.cast_nonneg (d N))]
    have hFmeas : Measurable (fun Z : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
        (Real.sqrt (d N))⁻¹ *
          (WithLp.ofLp (m.u N) ⬝ᵥ (Z *ᵥ WithLp.ofLp (m.v N)))) :=
      measurable_const.mul (measurable_dotProduct_mulVec _ _)
    have hHint := NoiseLaw.integrable_sq_dotProduct_mulVec hν
      (WithLp.ofLp (m.u N)) (WithLp.ofLp (m.v N))
    have hHbound := NoiseLaw.integral_sq_dotProduct_mulVec (ν := ν) hν hau hav
    have heq : (fun Z : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
        ((Real.sqrt (d N))⁻¹ * (WithLp.ofLp (m.u N) ⬝ᵥ (Z *ᵥ WithLp.ofLp (m.v N)))) ^ 2)
        = fun Z => ((Real.sqrt (d N))⁻¹) ^ 2
          * (WithLp.ofLp (m.u N) ⬝ᵥ (Z *ᵥ WithLp.ofLp (m.v N))) ^ 2 := by
      funext Z
      rw [mul_pow]
    have hFint : Integrable (fun Z : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
        ((Real.sqrt (d N))⁻¹ * (WithLp.ofLp (m.u N) ⬝ᵥ (Z *ᵥ WithLp.ofLp (m.v N)))) ^ 2)
        (noiseMatrix ν (n N) (d N)) := by
      rw [heq]
      exact hHint.const_mul _
    have hFbound : ∫ Z, ((Real.sqrt (d N))⁻¹
          * (WithLp.ofLp (m.u N) ⬝ᵥ (Z *ᵥ WithLp.ofLp (m.v N)))) ^ 2
        ∂(noiseMatrix ν (n N) (d N)) = 1 / (d N : ℝ) := by
      rw [heq, integral_const_mul, hHbound, mul_one, pow_two, hdd, one_div]
    have hFabsmeas := hFmeas.abs
    have hFabsint : Integrable (fun Z : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
        |(Real.sqrt (d N))⁻¹ * (WithLp.ofLp (m.u N) ⬝ᵥ (Z *ᵥ WithLp.ofLp (m.v N)))| ^ 2)
        (noiseMatrix ν (n N) (d N)) := by
      simpa only [sq_abs] using hFint
    have hFabsC : ∫ Z, |(Real.sqrt (d N))⁻¹
          * (WithLp.ofLp (m.u N) ⬝ᵥ (Z *ᵥ WithLp.ofLp (m.v N)))| ^ 2
        ∂(noiseMatrix ν (n N) (d N)) ≤ 1 / (d N : ℝ) := by
      simpa only [sq_abs] using hFbound.le
    have hset : MeasurableSet {Z : Matrix (Fin (n N)) (Fin (d N)) ℝ |
        ε ≤ |(Real.sqrt (d N))⁻¹
          * (WithLp.ofLp (m.u N) ⬝ᵥ (Z *ᵥ WithLp.ofLp (m.v N)))|} :=
      measurableSet_le measurable_const hFabsmeas
    have hchb := Cheb.meas_ge_le_of_integral_sq _ hFabsmeas hFabsint hε hFabsC
    have hmeq := (hG N).measure_eq hset
    have hpt : ∀ ω, (Real.sqrt (d N))⁻¹
          * (WithLp.ofLp (m.u N) ⬝ᵥ ((m.Z N ω) *ᵥ WithLp.ofLp (m.v N)))
        = WithLp.ofLp (m.v N) ⬝ᵥ m.gvec N ω := by
      intro ω
      rw [m.gvec_eq_smul N ω, dotProduct_smul, smul_eq_mul, Matrix.dotProduct_transpose_mulVec]
    have hseteq : {ω | ε ≤ |WithLp.ofLp (m.v N) ⬝ᵥ m.gvec N ω - 0|}
        = {ω | ε ≤ |(Real.sqrt (d N))⁻¹
            * (WithLp.ofLp (m.u N) ⬝ᵥ ((m.Z N ω) *ᵥ WithLp.ofLp (m.v N)))|} := by
      ext ω
      simp only [Set.mem_ofPred_eq, sub_zero, hpt]
    rw [hseteq, hmeq]
    exact hchb
  have h0 := Cheb.tendstoInProb_of_meas_le hK hbnd
  simpa using h0

/-- `q ⬝ q → θ² + 1`. General-law twin of `tendstoInProb_dotProduct_qvec` (`RMT/R6.lean:830`). -/
theorem tendstoInProb_dotProduct_qvec_general (m : SpikedModel μ n d) {ν : Measure ℝ}
    (hν : NoiseLaw ν) (hG : m.GeneralNoise ν) (hd : Tendsto d atTop atTop) :
    TendstoInProb μ (fun N ω => m.qvec N ω ⬝ᵥ m.qvec N ω) (m.θ ^ 2 + 1) := by
  have hfun : ∀ (N : ℕ) (ω : Ω N), m.qvec N ω ⬝ᵥ m.qvec N ω
      = m.θ ^ 2 + 2 * m.θ * (WithLp.ofLp (m.v N) ⬝ᵥ m.gvec N ω)
        + m.gvec N ω ⬝ᵥ m.gvec N ω := by
    intro N ω
    change (m.θ • WithLp.ofLp (m.v N) + m.gvec N ω)
        ⬝ᵥ (m.θ • WithLp.ofLp (m.v N) + m.gvec N ω) = _
    simp only [add_dotProduct, dotProduct_add, smul_dotProduct, dotProduct_smul, smul_eq_mul,
      m.dotProduct_v_self N,
      dotProduct_comm (m.gvec N ω) (WithLp.ofLp (m.v N))]
    ring
  have h := ((TendstoInProb.const μ (m.θ ^ 2)).add
    ((m.tendstoInProb_dotProduct_v_gvec_general hν hG hd).const_mul (2 * m.θ))).add
      (m.tendstoInProb_dotProduct_gvec_general hν hG hd)
  have hlim : m.θ ^ 2 + 2 * m.θ * 0 + 1 = m.θ ^ 2 + 1 := by ring
  rw [hlim] at h
  simpa only [hfun] using h

/-! ### the six resolvent forms at complex z -/

/-- The six isotropic resolvent forms of `ResolventLimits`, at every complex `z` with
`Im z > 0`, on `W₀`, `v` and `g`. Same six functions as the `hforms` block of
`resolventLimits_of_gaussian` (`RMT/Sup.lean`, step 4). Until F35 (2026-09-09) this
structure lived in `RMT/General/Statements.lean`, the file that held the 12 target
statements of the plan while they carried `sorry`; that file is retired. -/
structure ResolventFormsC (m : SpikedModel μ n d) (c : ℝ) : Prop where
  vvC : ∀ z : ℂ, 0 < z.im → TendstoInProb μ
    (fun N ω => ‖R4C.qformC (m.W0 N ω) z (WithLp.ofLp (m.v N)) - MP.mC c z‖) 0
  ggC : ∀ z : ℂ, 0 < z.im → TendstoInProb μ
    (fun N ω => ‖R4C.qformC (m.W0 N ω) z (m.gvec N ω) - MP.mC c z‖) 0
  vgC : ∀ z : ℂ, 0 < z.im → TendstoInProb μ
    (fun N ω => ‖R4C.cformC (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω)‖) 0
  vv2C : ∀ z : ℂ, 0 < z.im → TendstoInProb μ
    (fun N ω => ‖R4C.qform2C (m.W0 N ω) z (WithLp.ofLp (m.v N)) - MP.mCDeriv c z‖) 0
  gg2C : ∀ z : ℂ, 0 < z.im → TendstoInProb μ
    (fun N ω => ‖R4C.qform2C (m.W0 N ω) z (m.gvec N ω) - MP.mCDeriv c z‖) 0
  vg2C : ∀ z : ℂ, 0 < z.im → TendstoInProb μ
    (fun N ω => ‖R4C.cform2C (m.W0 N ω) z (WithLp.ofLp (m.v N)) (m.gvec N ω)‖) 0

end SpikedModel
end StackedSVD
