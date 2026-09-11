/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Edge.Arith
import StackedSVD.RMT.General.Edge.Trunc
import StackedSVD.RMT.General.Edge.Trace
import StackedSVD.RMT.General.Edge.Compare
import StackedSVD.RMT.General.Edge.Excess
import StackedSVD.RMT.General.Edge.Gaussian
import StackedSVD.RMT.R0
import StackedSVD.RMT.R3
import StackedSVD.Prob.TendstoInProb

/-! # Stage 3, unit K (second half) and unit M: the Markov step

The two spectral contractions (`lamMax(W)^k ≤ trace (W^k)`, `‖Y‖² ≤ lamMax (YᵀY)`), the
rank-one downdate `‖E⊥‖ ≤ ‖E‖`, the bundled moment bound of unit M (units C, G, X and T at
the instantiation of the assembly, proved here), and the Markov step K8: with probability
tending to 1, `‖Ẑ‖²/d ≤ bulkEdge c + 2ε` (`notes/stage3_edge.md`, route C;
`notes/STAGE3_CAMPAIGN.md`). -/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix Matrix.Norms.L2Operator ENNReal NNReal

namespace StackedSVD
namespace Edge

/-- The Gram quadratic form: `x ⬝ᵥ (YᵀY x) = ‖Y x‖²`. -/
private theorem dotProduct_gram {n D : ℕ} (Y : Matrix (Fin n) (Fin D) ℝ) (x : Fin D → ℝ) :
    x ⬝ᵥ ((Yᵀ * Y) *ᵥ x) = (Y *ᵥ x) ⬝ᵥ (Y *ᵥ x) := by
  rw [← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec, Matrix.vecMul_transpose]

/-- A unit eigenvector at the top eigenvalue, in the eigenbasis of `R4.eigU`. -/
private theorem exists_top_eigvec {D : ℕ} (hD : 0 < D) {W : Matrix (Fin D) (Fin D) ℝ}
    (hW : W.IsHermitian) :
    ∃ y : Fin D → ℝ, y ⬝ᵥ y = 1 ∧ W *ᵥ y = lamMax W hW • y := by
  classical
  obtain ⟨j, hj⟩ := R4.exists_eigenvalues_eq_lamMax hW hD
  refine ⟨R4.eigU hW *ᵥ Pi.single j 1, ?_, ?_⟩
  · have hcoord : (R4.eigU hW)ᵀ *ᵥ (R4.eigU hW *ᵥ Pi.single j (1 : ℝ)) = Pi.single j 1 := by
      rw [Matrix.mulVec_mulVec, R4.transpose_eigU_mul, Matrix.one_mulVec]
    rw [← R4.dotProduct_transpose_eigU hW, hcoord]
    simp [dotProduct, Pi.single_apply]
  · have hcoord : (R4.eigU hW)ᵀ *ᵥ (R4.eigU hW *ᵥ Pi.single j (1 : ℝ)) = Pi.single j 1 := by
      rw [Matrix.mulVec_mulVec, R4.transpose_eigU_mul, Matrix.one_mulVec]
    have hconj : W *ᵥ (R4.eigU hW *ᵥ Pi.single j (1 : ℝ))
        = (R4.eigU hW * Matrix.diagonal hW.eigenvalues * (R4.eigU hW)ᵀ)
            *ᵥ (R4.eigU hW *ᵥ Pi.single j (1 : ℝ)) := by
      rw [R4.eigU_conj hW]
    have hfun : (fun a => hW.eigenvalues a * (Pi.single j (1 : ℝ) : Fin D → ℝ) a)
        = hW.eigenvalues j • (Pi.single j (1 : ℝ) : Fin D → ℝ) := by
      funext b
      by_cases hb : b = j
      · subst hb; simp
      · simp [hb]
    rw [hconj, R4.conj_mulVec, hcoord, hfun, Matrix.mulVec_smul, hj]

/-- `W^k y = λ^k y` at an eigenvector `y` of eigenvalue `λ`. -/
private theorem pow_mulVec_of_eigvec {D : ℕ} {W : Matrix (Fin D) (Fin D) ℝ} {y : Fin D → ℝ}
    {lam : ℝ} (h : W *ᵥ y = lam • y) (k : ℕ) : (W ^ k) *ᵥ y = lam ^ k • y := by
  induction k with
  | zero => simp
  | succ p ih =>
      rw [pow_succ, ← Matrix.mulVec_mulVec, h, Matrix.mulVec_smul, ih, smul_smul, pow_succ]
      ring_nf

/-- K5. `lamMax(W)^k ≤ trace (W^k)` for a positive semidefinite `W`. Through the spectral
theorem: `W = U D U⋆`, `trace (W^k) = ∑ λ_i^k`, every term nonnegative. (`Spectral.lean` has
a private `quadForm_pow`; it stays private, decision D-M of the campaign note.) -/
theorem lamMax_pow_le_trace {D : ℕ} (hD : 0 < D) (W : Matrix (Fin D) (Fin D) ℝ)
    (hW : W.PosSemidef) (k : ℕ) :
    lamMax W hW.isHermitian ^ k ≤ Matrix.trace (W ^ k) := by
  classical
  obtain ⟨y, hyy, hy⟩ := exists_top_eigvec hD hW.isHermitian
  have hWk : (W ^ k).PosSemidef := hW.pow k
  have h1 : y ⬝ᵥ ((W ^ k) *ᵥ y) ≤ lamMax (W ^ k) hWk.isHermitian * (y ⬝ᵥ y) :=
    R4.dotProduct_mulVec_le_lamMax hWk.isHermitian y
  have h2 : y ⬝ᵥ ((W ^ k) *ᵥ y) = lamMax W hW.isHermitian ^ k := by
    rw [pow_mulVec_of_eigvec hy k, dotProduct_smul, smul_eq_mul, hyy, mul_one]
  rw [hyy, mul_one, h2] at h1
  obtain ⟨i0, hi0⟩ := R4.exists_eigenvalues_eq_lamMax hWk.isHermitian hD
  have h3 : Matrix.trace (W ^ k) = ∑ i, hWk.isHermitian.eigenvalues i := by
    simpa using hWk.isHermitian.trace_eq_sum_eigenvalues
  have h4 : lamMax (W ^ k) hWk.isHermitian ≤ Matrix.trace (W ^ k) := by
    rw [h3, ← hi0]
    exact Finset.single_le_sum (fun i _ => hWk.eigenvalues_nonneg i) (Finset.mem_univ i0)
  linarith

/-- K6. The operator norm of the noise against the top eigenvalue of the Gram matrix:
`‖Y‖² = lamMax (YᵀY)`, in the direction the assembly needs. -/
theorem opNorm_sq_le_lamMax_gram {n D : ℕ} (hD : 0 < D) (Y : Matrix (Fin n) (Fin D) ℝ) :
    ‖Y‖ ^ 2 ≤ lamMax (Yᵀ * Y) (isHermitian_transpose_mul_self Y) := by
  have hray : ∀ x : Fin D → ℝ, (Y *ᵥ x) ⬝ᵥ (Y *ᵥ x)
      ≤ lamMax (Yᵀ * Y) (isHermitian_transpose_mul_self Y) * (x ⬝ᵥ x) := by
    intro x
    rw [← dotProduct_gram Y x]
    exact R4.dotProduct_mulVec_le_lamMax _ x
  have hL0 : (0 : ℝ) ≤ lamMax (Yᵀ * Y) (isHermitian_transpose_mul_self Y) := by
    have h := hray (Pi.single (⟨0, hD⟩ : Fin D) 1)
    have hxx : (Pi.single (⟨0, hD⟩ : Fin D) (1 : ℝ)) ⬝ᵥ (Pi.single (⟨0, hD⟩ : Fin D) 1) = 1 := by
      simp [dotProduct, Pi.single_apply]
    have hnn : (0 : ℝ) ≤ (Y *ᵥ Pi.single (⟨0, hD⟩ : Fin D) 1) ⬝ᵥ
        (Y *ᵥ Pi.single (⟨0, hD⟩ : Fin D) 1) :=
      Finset.sum_nonneg fun i _ => mul_self_nonneg _
    rw [hxx, mul_one] at h
    linarith
  have hnorm : ‖Y‖ ≤ Real.sqrt (lamMax (Yᵀ * Y) (isHermitian_transpose_mul_self Y)) := by
    rw [Matrix.l2_opNorm_def]
    refine ContinuousLinearMap.opNorm_le_bound _ (Real.sqrt_nonneg _) fun x => ?_
    have happ : ((Matrix.toEuclideanLin.trans LinearMap.toContinuousLinearMap) Y) x
        = WithLp.toLp 2 (Y *ᵥ WithLp.ofLp x) := rfl
    rw [happ]
    have hA : ‖(WithLp.toLp 2 (Y *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin n))‖ ^ 2
        = (Y *ᵥ WithLp.ofLp x) ⬝ᵥ (Y *ᵥ WithLp.ofLp x) := by
      rw [EuclideanSpace.real_norm_sq_eq]
      simp [dotProduct, sq]
    have hB : ‖x‖ ^ 2 = (WithLp.ofLp x) ⬝ᵥ (WithLp.ofLp x) := by
      rw [EuclideanSpace.real_norm_sq_eq]
      simp [dotProduct, sq]
    have hsq : ‖(WithLp.toLp 2 (Y *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin n))‖ ^ 2
        ≤ lamMax (Yᵀ * Y) (isHermitian_transpose_mul_self Y) * ‖x‖ ^ 2 := by
      rw [hA, hB]
      exact hray _
    calc ‖(WithLp.toLp 2 (Y *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin n))‖
        = Real.sqrt (‖(WithLp.toLp 2 (Y *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin n))‖ ^ 2) :=
          (Real.sqrt_sq (norm_nonneg _)).symm
      _ ≤ Real.sqrt (lamMax (Yᵀ * Y) (isHermitian_transpose_mul_self Y) * ‖x‖ ^ 2) :=
          Real.sqrt_le_sqrt hsq
      _ = Real.sqrt (lamMax (Yᵀ * Y) (isHermitian_transpose_mul_self Y)) * ‖x‖ := by
          rw [Real.sqrt_mul hL0, Real.sqrt_sq (norm_nonneg _)]
  have hfin := pow_le_pow_left₀ (norm_nonneg Y) hnorm 2
  rwa [Real.sq_sqrt hL0] at hfin

/-- K7. The rank-one downdate lowers the operator norm: `‖E⊥‖ ≤ ‖E‖`, since
`E⊥ = (I - u uᵀ) E` (choice 4 of the note). -/
theorem opNorm_Eperp_le {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
    {n d : ℕ → ℕ} (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) :
    ‖m.Eperp N ω‖ ≤ ‖m.E N ω‖ := by
  have huu : WithLp.ofLp (m.u N) ⬝ᵥ WithLp.ofLp (m.u N) = 1 := m.dotProduct_u_self N
  -- the rank-one downdate acts as the orthogonal projector `1 - u uᵀ` on every image vector
  have hmv : ∀ v : Fin (d N) → ℝ, (m.Eperp N ω) *ᵥ v
      = (m.E N ω *ᵥ v) - (WithLp.ofLp (m.u N) ⬝ᵥ (m.E N ω *ᵥ v)) • WithLp.ofLp (m.u N) := by
    intro v
    have hEp : m.Eperp N ω
        = m.E N ω - Matrix.vecMulVec (WithLp.ofLp (m.u N)) ((m.E N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)) :=
      rfl
    rw [hEp, Matrix.sub_mulVec]
    congr 1
    funext i
    have hdot : ((m.E N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)) ⬝ᵥ v
        = WithLp.ofLp (m.u N) ⬝ᵥ (m.E N ω *ᵥ v) := by
      rw [Matrix.mulVec_transpose, Matrix.dotProduct_mulVec]
    calc (Matrix.vecMulVec (WithLp.ofLp (m.u N)) ((m.E N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)) *ᵥ v) i
        = WithLp.ofLp (m.u N) i * (((m.E N ω)ᵀ *ᵥ WithLp.ofLp (m.u N)) ⬝ᵥ v) := by
          simp only [Matrix.mulVec, dotProduct, Matrix.vecMulVec_apply, Matrix.transpose_apply,
            Finset.mul_sum]
          refine Finset.sum_congr rfl fun x _ => ?_
          rw [← Finset.mul_sum, mul_assoc]
      _ = ((WithLp.ofLp (m.u N) ⬝ᵥ (m.E N ω *ᵥ v)) • WithLp.ofLp (m.u N)) i := by
          rw [hdot]; simp [mul_comm]
  have hkey : ∀ v : Fin (d N) → ℝ,
      ((m.Eperp N ω) *ᵥ v) ⬝ᵥ ((m.Eperp N ω) *ᵥ v) ≤ (m.E N ω *ᵥ v) ⬝ᵥ (m.E N ω *ᵥ v) := by
    intro v
    rw [hmv v]
    have hexp : ((m.E N ω *ᵥ v) - (WithLp.ofLp (m.u N) ⬝ᵥ (m.E N ω *ᵥ v)) • WithLp.ofLp (m.u N))
          ⬝ᵥ ((m.E N ω *ᵥ v)
            - (WithLp.ofLp (m.u N) ⬝ᵥ (m.E N ω *ᵥ v)) • WithLp.ofLp (m.u N))
        = (m.E N ω *ᵥ v) ⬝ᵥ (m.E N ω *ᵥ v)
          - (WithLp.ofLp (m.u N) ⬝ᵥ (m.E N ω *ᵥ v)) ^ 2 := by
      simp only [sub_dotProduct, dotProduct_sub, smul_dotProduct, dotProduct_smul,
        smul_eq_mul, dotProduct_comm (m.E N ω *ᵥ v) (WithLp.ofLp (m.u N)), huu]
      ring
    rw [hexp]
    nlinarith [sq_nonneg (WithLp.ofLp (m.u N) ⬝ᵥ (m.E N ω *ᵥ v))]
  rw [Matrix.l2_opNorm_def]
  refine ContinuousLinearMap.opNorm_le_bound _ (norm_nonneg _) fun x => ?_
  have happ : ((Matrix.toEuclideanLin.trans LinearMap.toContinuousLinearMap) (m.Eperp N ω)) x
      = WithLp.toLp 2 (m.Eperp N ω *ᵥ WithLp.ofLp x) := rfl
  rw [happ]
  have hsq : ∀ {p : ℕ} (z : Fin p → ℝ),
      ‖(WithLp.toLp 2 z : EuclideanSpace ℝ (Fin p))‖ ^ 2 = z ⬝ᵥ z := by
    intro p z
    rw [EuclideanSpace.real_norm_sq_eq]
    simp [dotProduct, sq]
  have hle : ‖(WithLp.toLp 2 (m.Eperp N ω *ᵥ WithLp.ofLp x) :
        EuclideanSpace ℝ (Fin (n N)))‖
      ≤ ‖(WithLp.toLp 2 (m.E N ω *ᵥ WithLp.ofLp x) : EuclideanSpace ℝ (Fin (n N)))‖ := by
    have h := hkey (WithLp.ofLp x)
    rw [← hsq (m.Eperp N ω *ᵥ WithLp.ofLp x), ← hsq (m.E N ω *ᵥ WithLp.ofLp x)] at h
    have h2 := Real.sqrt_le_sqrt h
    rwa [Real.sqrt_sq (norm_nonneg _), Real.sqrt_sq (norm_nonneg _)] at h2
  refine hle.trans ?_
  exact ((Matrix.toEuclideanLin.trans LinearMap.toContinuousLinearMap) (m.E N ω)).le_opNorm x

/-- M. `∫ trace ((ẐᵀẐ)^k) ≤ 4 d (d (bulkEdge c + ε))^k` at `k = k_N`, for `N` large.
Units C, G, X and T at the instantiation of the assembly. -/
theorem trunc_trace_le {c : ℝ} (hc : 0 < c) {ν : Measure ℝ} (hν : NoiseLaw ν) {n d : ℕ → ℕ}
    (hn : ∀ N, 0 < n N) (hd : ∀ N, 0 < d N) (hdtop : Tendsto d atTop atTop)
    (hcN : Tendsto (fun N => (n N : ℝ) / d N) atTop (𝓝 c))
    {a : ℝ} (ha : 0 < a) (ha2 : a ≤ 1 / 2) {C : ℝ} (hC : 0 < C) {ε : ℝ} (hε : 0 < ε) :
    ∀ᶠ N in atTop,
      ∫ Y, Matrix.trace ((Yᵀ * Y) ^ (momOrder C (d N)))
          ∂(noiseMatrix (truncLaw ν (truncLevel a (d N))) (n N) (d N))
        ≤ 4 * (d N : ℝ) * ((d N : ℝ) * (bulkEdge c + ε)) ^ (momOrder C (d N)) := by
  filter_upwards [gaussian_trace_le hc hn hd hdtop hcN hC hε,
    eventually_fkFactor_le hc hdtop hcN ha hC,
    eventually_one_le_momOrder hdtop hC] with N h3 hf hk1
  set k : ℕ := momOrder C (d N) with hkdef
  set T : ℝ := truncLevel a (d N) with hTdef
  have hT : 0 < T := truncLevel_pos (hd N)
  -- `1 ≤ T = d^(1/2 - a)` since `1 ≤ d` and `a ≤ 1/2`; unit X needs `1 ≤ K = 2T`
  have hT1 : (1 : ℝ) ≤ T := by
    rw [hTdef, truncLevel]
    exact Real.one_le_rpow (Nat.one_le_cast.mpr (hd N)) (by linarith)
  have hK0 : (0 : ℝ) ≤ 2 * T := by linarith
  have hK : (1 : ℝ) ≤ 2 * T := by linarith
  have hrho : TruncNoiseLaw (truncLaw ν T) (2 * T) := truncNoiseLaw_truncLaw hν hT
  have hGnn : (0 : ℝ) ≤ ∫ Y, Matrix.trace ((Yᵀ * Y) ^ k) ∂(gaussianMatrix (n N) (d N)) :=
    integral_trace_gaussian_nonneg _ _ _
  have h1 := trace_le_gaussian_add_excess hK0 hrho (n := n N) (d := d N) (k := k) hk1
  have h2 := excessSum_le hK hrho hk1 (hn N) (hd N) hf
  have hfnn : (0 : ℝ) ≤ fkFactor (2 * T) (n N) (d N) k := by
    rw [fkFactor]; positivity
  nlinarith [h1, h2, h3, hGnn, hfnn, hf]

/-- K8, the Markov step. Every hypothesis is the conclusion of an earlier unit, at the
instantiation the assembly uses. -/
theorem tendsto_measure_trunc_opNorm {c : ℝ} (hc : 0 < c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    {n d : ℕ → ℕ} (hn : ∀ N, 0 < n N) (hd : ∀ N, 0 < d N)
    (hdtop : Tendsto d atTop atTop) (hcN : Tendsto (fun N => (n N : ℝ) / d N) atTop (𝓝 c))
    {a : ℝ} (ha : 0 < a) (ha2 : a ≤ 1 / 2) {ε : ℝ} (hε : 0 < ε) :
    Tendsto (fun N => (noiseMatrix ν (n N) (d N))
        {Z | ‖truncMat ν (truncLevel a (d N)) Z‖ ^ 2 / (d N : ℝ) ≤ bulkEdge c + 2 * ε})
      atTop (𝓝 1) := by
  classical
  have : IsProbabilityMeasure ν := hν.prob
  have hbe : (0 : ℝ) < bulkEdge c := by rw [bulkEdge]; positivity
  have hb1 : (0 : ℝ) < bulkEdge c + ε := by linarith
  have hb2 : (0 : ℝ) < bulkEdge c + 2 * ε := by linarith
  have hCpos : 0 < markovConst c ε := markovConst_pos hc hε
  -- the bad trace event, on the truncated law
  obtain ⟨S, hS⟩ : ∃ S : (N : ℕ) → Set (Matrix (Fin (n N)) (Fin (d N)) ℝ), ∀ N, S N =
      {Y | ((d N : ℝ) * (bulkEdge c + 2 * ε)) ^ (momOrder (markovConst c ε) (d N))
        ≤ Matrix.trace ((Yᵀ * Y) ^ (momOrder (markovConst c ε) (d N)))} :=
    ⟨_, fun _ => rfl⟩
  have hSmeas : ∀ N, MeasurableSet (S N) := by
    intro N
    rw [hS N]
    exact measurableSet_le measurable_const
      (measurable_trace_gram_pow (n N) (d N) (momOrder (markovConst c ε) (d N)))
  refine tendsto_measure_one_of_bad
    (Ω := fun N => Matrix (Fin (n N)) (Fin (d N)) ℝ)
    (μ := fun N => noiseMatrix ν (n N) (d N))
    (t := fun N => {Z | ‖truncMat ν (truncLevel a (d N)) Z‖ ^ 2 / (d N : ℝ)
      ≤ bulkEdge c + 2 * ε})
    (s := fun N => truncMat ν (truncLevel a (d N)) ⁻¹' (S N)) ?_ ?_
  · -- the two spectral contractions put the bad norm event inside the bad trace event
    intro N Z hZ
    have hd0 : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast hd N
    have hZ' : bulkEdge c + 2 * ε
        < ‖truncMat ν (truncLevel a (d N)) Z‖ ^ 2 / (d N : ℝ) := by
      simpa using hZ
    have h1 : (d N : ℝ) * (bulkEdge c + 2 * ε)
        ≤ ‖truncMat ν (truncLevel a (d N)) Z‖ ^ 2 := by
      rw [lt_div_iff₀ hd0] at hZ'
      linarith
    have hpsd : ((truncMat ν (truncLevel a (d N)) Z)ᵀ
        * truncMat ν (truncLevel a (d N)) Z).PosSemidef := by
      have h0 := Matrix.posSemidef_conjTranspose_mul_self (truncMat ν (truncLevel a (d N)) Z)
      rwa [Matrix.conjTranspose_eq_transpose_of_trivial] at h0
    have h2 := opNorm_sq_le_lamMax_gram (hd N) (truncMat ν (truncLevel a (d N)) Z)
    have h3 : ((d N : ℝ) * (bulkEdge c + 2 * ε)) ^ (momOrder (markovConst c ε) (d N))
        ≤ (lamMax ((truncMat ν (truncLevel a (d N)) Z)ᵀ * truncMat ν (truncLevel a (d N)) Z)
            (isHermitian_transpose_mul_self _)) ^ (momOrder (markovConst c ε) (d N)) :=
      pow_le_pow_left₀ (by positivity) (by linarith) _
    have h4 := lamMax_pow_le_trace (hd N)
      ((truncMat ν (truncLevel a (d N)) Z)ᵀ * truncMat ν (truncLevel a (d N)) Z) hpsd
      (momOrder (markovConst c ε) (d N))
    change truncMat ν (truncLevel a (d N)) Z ∈ S N
    rw [hS N]
    exact h3.trans h4
  · -- Markov, unit M, and K4
    have htrans : ∀ N, (noiseMatrix ν (n N) (d N)) (truncMat ν (truncLevel a (d N)) ⁻¹' (S N))
        = (noiseMatrix (truncLaw ν (truncLevel a (d N))) (n N) (d N)) (S N) := fun N =>
      measure_truncMat_preimage ν (truncLevel a (d N)) (n N) (d N) (hSmeas N)
    simp only [htrans]
    have hreal : Tendsto (fun N =>
        (noiseMatrix (truncLaw ν (truncLevel a (d N))) (n N) (d N)).real (S N)) atTop (𝓝 0) := by
      refine squeeze_zero' (Eventually.of_forall fun N => measureReal_nonneg) ?_
        (tendsto_markov_zero hc hε hdtop)
      filter_upwards [trunc_trace_le hc hν hn hd hdtop hcN ha ha2 hCpos hε] with N hMN
      have hd0 : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast hd N
      have hdne : (d N : ℝ) ≠ 0 := ne_of_gt hd0
      have hT : 0 < truncLevel a (d N) := truncLevel_pos (hd N)
      have hrho : TruncNoiseLaw (truncLaw ν (truncLevel a (d N)))
          (2 * truncLevel a (d N)) := truncNoiseLaw_truncLaw hν hT
      have hint : Integrable (fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
          Matrix.trace ((Yᵀ * Y) ^ (momOrder (markovConst c ε) (d N))))
          (noiseMatrix (truncLaw ν (truncLevel a (d N))) (n N) (d N)) :=
        integrable_trace_gram_pow (by linarith) hrho
      have hmk : ((d N : ℝ) * (bulkEdge c + 2 * ε)) ^ (momOrder (markovConst c ε) (d N))
          * (noiseMatrix (truncLaw ν (truncLevel a (d N))) (n N) (d N)).real (S N)
          ≤ ∫ Y, Matrix.trace ((Yᵀ * Y) ^ (momOrder (markovConst c ε) (d N)))
              ∂(noiseMatrix (truncLaw ν (truncLevel a (d N))) (n N) (d N)) := by
        rw [hS N]
        exact mul_meas_ge_le_integral_of_nonneg
          (ae_of_all _ fun Y => trace_gram_pow_nonneg Y) hint _
      have hchain := hmk.trans hMN
      have hthr : (0 : ℝ) < ((d N : ℝ) * (bulkEdge c + 2 * ε))
          ^ (momOrder (markovConst c ε) (d N)) := pow_pos (mul_pos hd0 hb2) _
      have hcomm : (noiseMatrix (truncLaw ν (truncLevel a (d N))) (n N) (d N)).real (S N)
          * ((d N : ℝ) * (bulkEdge c + 2 * ε)) ^ (momOrder (markovConst c ε) (d N))
          ≤ 4 * (d N : ℝ)
              * ((d N : ℝ) * (bulkEdge c + ε)) ^ (momOrder (markovConst c ε) (d N)) := by
        rw [mul_comm]; exact hchain
      have hdiv := (le_div_iff₀ hthr).mpr hcomm
      have hcanc : 4 * (d N : ℝ)
            * ((d N : ℝ) * (bulkEdge c + ε)) ^ (momOrder (markovConst c ε) (d N))
            / ((d N : ℝ) * (bulkEdge c + 2 * ε)) ^ (momOrder (markovConst c ε) (d N))
          = 4 * (d N : ℝ)
            * ((bulkEdge c + ε) / (bulkEdge c + 2 * ε)) ^ (momOrder (markovConst c ε) (d N)) := by
        have h1 : ((d N : ℝ)) ^ (momOrder (markovConst c ε) (d N)) ≠ 0 := pow_ne_zero _ hdne
        have h2 : (bulkEdge c + 2 * ε) ^ (momOrder (markovConst c ε) (d N)) ≠ 0 :=
          pow_ne_zero _ (ne_of_gt hb2)
        rw [mul_pow, mul_pow, div_pow]
        field_simp
      rwa [hcanc] at hdiv
    have hgoal := (ENNReal.continuous_ofReal.tendsto 0).comp hreal
    simp only [Function.comp_def, ENNReal.ofReal_zero] at hgoal
    refine hgoal.congr fun N => ?_
    rw [measureReal_def, ENNReal.ofReal_toReal (measure_ne_top _ _)]

end Edge
end StackedSVD
