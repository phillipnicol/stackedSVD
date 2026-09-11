/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.Subspace

/-!
# The shared stack interface `RankRStack`

Task A1 of `notes/archive/rankr_plan_A.md` section 2 (the `[AUDIT]` paragraph), answering open
question 1 of `notes/archive/rankr_plan_B.md`. The Gaussian discharge chain of the stacked table
(`exists_rankR_split`, `resolventLimitsR_of_gaussian`,
`subspaceLaw_of_gaussian_supercritical`) reads six objects of `UnalignedModel`:
`signalPart`, `spikeMat`, `coreEig`, `stackGram`, `stackE`, `spikeVec`. None of the content
depends on `r_i = 1`. This file states the data those objects are built from, so that a
general-`r_i` stack (`UnalignedModelR` with `stackXG` / `coreG` / `spikeVecG`, plan B
section 4) instantiates the same chain.

The data are: a shared orthonormal `V ∈ ℝ^{d × r}`, a signal factor `A ∈ ℝ^{ns × r}` whose
Gram is an `N`-free core matrix `C`, an unscaled noise `Zu`, its scaling `E = d^{-1/2} Zu`,
and the table `X = A Vᵀ + E`. `E` and `X` are fields, not definitions, so that a model whose
own `stackE` and `stackX` are built by another route (a block reindexing, say) instantiates
the structure with `rfl` bridges and keeps every derived object definitionally equal.

Nothing here is random: `GaussianNoise` is a predicate, as `SpikedModel.GaussianNoise` is,
and every theorem of the chain takes it as a hypothesis.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

/-- **The shared stack data of the rank-`r` chain.** One table `X = A Vᵀ + E` per `N`, with
`E = d^{-1/2} Zu`, orthonormal `V`, and an `N`-free core `C = AᵀA`.

`UnalignedModel` (`r_i = 1`) instantiates it through `UnalignedModel.toStack`; a general-`r_i`
model instantiates it with the same fields and a different signal factor. -/
structure RankRStack {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    (μ : ∀ N, Measure (Ω N)) (ns : ℕ → ℕ) (d : ℕ → ℕ) (r : ℕ) where
  /-- the shared ambient subspace, as a matrix with orthonormal columns -/
  V : (N : ℕ) → Matrix (Fin (d N)) (Fin r) ℝ
  /-- the signal factor `A` of the stack, `E[X] = A Vᵀ` -/
  A : (N : ℕ) → Matrix (Fin (ns N)) (Fin r) ℝ
  /-- the core matrix `C = AᵀA`, the same at every `N` -/
  core : Matrix (Fin r) (Fin r) ℝ
  /-- the unscaled noise, Gaussian under `GaussianNoise` -/
  Zu : (N : ℕ) → Ω N → Matrix (Fin (ns N)) (Fin (d N)) ℝ
  /-- the scaled noise `E = d^{-1/2} Zu` -/
  E : (N : ℕ) → Ω N → Matrix (Fin (ns N)) (Fin (d N)) ℝ
  /-- the table itself, `X = A Vᵀ + E` -/
  X : (N : ℕ) → Ω N → Matrix (Fin (ns N)) (Fin (d N)) ℝ
  hV : ∀ N, (V N)ᵀ * V N = 1
  hcore : ∀ N, (A N)ᵀ * A N = core
  hE : ∀ N ω, E N ω = (Real.sqrt (d N))⁻¹ • Zu N ω
  hX : ∀ N ω, X N ω = A N * (V N)ᵀ + E N ω
  hd : ∀ N, 0 < d N
  hZmeas : ∀ N, Measurable (Zu N)

namespace RankRStack

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-- `assum:general_noise` for the stack: the unscaled noise is a canonical Gaussian matrix.
On `UnalignedModel` this is `hasLaw_stackZu`, that is `MultiTableModel.stack_law`. -/
def GaussianNoise (s : RankRStack μ ns d r) : Prop :=
  ∀ N, HasLaw (s.Zu N) (gaussianMatrix (ns N) (d N)) (μ N)

/-! ### 1. The table and its Gram matrix -/

/-- The mean of the table, `E[X] = A Vᵀ`. -/
noncomputable def signalPart (s : RankRStack μ ns d r) (N : ℕ) :
    Matrix (Fin (ns N)) (Fin (d N)) ℝ := s.A N * (s.V N)ᵀ

/-- Gram matrix of the table. -/
noncomputable def gram (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (d N)) (Fin (d N)) ℝ := (s.X N ω)ᵀ * s.X N ω

theorem isHermitian_gram (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N) :
    (s.gram N ω).IsHermitian :=
  isHermitian_transpose_mul_self (s.X N ω)

theorem X_eq (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N) :
    s.X N ω = s.signalPart N + s.E N ω := s.hX N ω

/-! ### 2. The core matrix and its spectrum -/

theorem isHermitian_core (s : RankRStack μ ns d r) : s.core.IsHermitian := by
  rw [← s.hcore 0]
  exact isHermitian_transpose_mul_self _

/-- `λ_j(C)`, at Mathlib's unsorted index, as `UnalignedModel.coreEig` is. -/
noncomputable def coreEig (s : RankRStack μ ns d r) (j : Fin r) : ℝ :=
  s.isHermitian_core.eigenvalues j

/-- `0 ≤ λ_j(C)`: the core matrix is a Gram matrix. -/
theorem coreEig_nonneg (s : RankRStack μ ns d r) (j : Fin r) : 0 ≤ s.coreEig j := by
  have hpsd : s.core.PosSemidef := by
    rw [← s.hcore 0]
    simpa using Matrix.posSemidef_conjTranspose_mul_self (s.A 0)
  exact hpsd.eigenvalues_nonneg j

/-- `Q_C`, the matrix whose columns are the eigenvectors of `C`. -/
noncomputable def coreEigMat (s : RankRStack μ ns d r) : Matrix (Fin r) (Fin r) ℝ :=
  (s.isHermitian_core.eigenvectorUnitary : Matrix (Fin r) (Fin r) ℝ)

theorem coreEigMat_apply (s : RankRStack μ ns d r) (k j : Fin r) :
    s.coreEigMat k j = s.isHermitian_core.eigenvectorBasis j k := rfl

theorem coreEigMat_transpose_mul_self (s : RankRStack μ ns d r) :
    (s.coreEigMat)ᵀ * s.coreEigMat = 1 := by
  have h : star (s.isHermitian_core.eigenvectorUnitary : Matrix (Fin r) (Fin r) ℝ) *
      (s.isHermitian_core.eigenvectorUnitary : Matrix (Fin r) (Fin r) ℝ) = 1 :=
    Unitary.coe_star_mul_self _
  have hst : (s.coreEigMat)ᵀ
      = star (s.isHermitian_core.eigenvectorUnitary : Matrix (Fin r) (Fin r) ℝ) := by
    ext k l
    simp [coreEigMat, Matrix.star_apply]
  rw [hst]
  exact h

theorem coreEigMat_mul_transpose_self (s : RankRStack μ ns d r) :
    s.coreEigMat * (s.coreEigMat)ᵀ = 1 :=
  mul_eq_one_comm.1 s.coreEigMat_transpose_mul_self

/-- `C = Q_C Λ Q_Cᵀ`, the spectral theorem in the matrix form the factorization takes. -/
theorem core_eq_conj (s : RankRStack μ ns d r) :
    s.core = s.coreEigMat * Matrix.diagonal s.coreEig * (s.coreEigMat)ᵀ := by
  ext k l
  conv_lhs => rw [s.isHermitian_core.spectral_theorem]
  simp only [Unitary.conjStarAlgAut_apply, Matrix.mul_apply, Matrix.diagonal_apply,
    Matrix.star_apply, Matrix.transpose_apply, Matrix.IsHermitian.eigenvectorUnitary_apply,
    star_trivial, Function.comp_apply, RCLike.ofReal_real_eq_id, id_eq, coreEigMat, coreEig,
    mul_ite, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, if_true]

/-! ### 3. The spike directions -/

/-- `V Q_C`, the `d × r` matrix whose columns are the spike directions `V q_j`. -/
noncomputable def spikeMat (s : RankRStack μ ns d r) (N : ℕ) :
    Matrix (Fin (d N)) (Fin r) ℝ := s.V N * s.coreEigMat

/-- `V q_j`, the `j`-th population right singular vector of the table. -/
noncomputable def spikeVec (s : RankRStack μ ns d r) (j : Fin r) (N : ℕ) :
    EuclideanSpace ℝ (Fin (d N)) :=
  WithLp.toLp 2 ((s.V N).mulVec (WithLp.ofLp (s.isHermitian_core.eigenvectorBasis j)))

theorem spikeMat_apply (s : RankRStack μ ns d r) (N : ℕ) (k : Fin (d N)) (j : Fin r) :
    s.spikeMat N k j = WithLp.ofLp (s.spikeVec j N) k := by
  simp [spikeMat, coreEigMat, spikeVec, Matrix.mul_apply, Matrix.mulVec, dotProduct]

/-- The spike directions are orthonormal: `(V Q_C)ᵀ (V Q_C) = 1`. -/
theorem spikeMat_transpose_mul_self (s : RankRStack μ ns d r) (N : ℕ) :
    (s.spikeMat N)ᵀ * s.spikeMat N = 1 := by
  rw [spikeMat, Matrix.transpose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc (s.V N)ᵀ,
    s.hV N, Matrix.one_mul, s.coreEigMat_transpose_mul_self]

/-- The spike directions are orthonormal, in the inner-product form the interface takes. -/
theorem inner_spikeVec (s : RankRStack μ ns d r) (N : ℕ) (a b : Fin r) :
    ⟪s.spikeVec a N, s.spikeVec b N⟫_ℝ = if a = b then 1 else 0 := by
  have h := congrFun (congrFun (s.spikeMat_transpose_mul_self N) a) b
  rw [Matrix.mul_apply] at h
  simp only [Matrix.transpose_apply, Matrix.one_apply] at h
  rw [inner_euclidean_eq_dotProduct]
  calc WithLp.ofLp (s.spikeVec a N) ⬝ᵥ WithLp.ofLp (s.spikeVec b N)
      = ∑ l, s.spikeMat N l a * s.spikeMat N l b := by
        simp only [dotProduct, s.spikeMat_apply N]
    _ = if a = b then 1 else 0 := h

/-! ### 4. The objects of the split -/

/-- `G = Eᵀ U` on the unscaled noise: `G = d^{-1/2} (Uᵀ Zu)ᵀ`. -/
theorem E_transpose_mul (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (ns N)) (Fin r) ℝ) :
    (s.E N ω)ᵀ * U = (Real.sqrt (d N))⁻¹ • ((Uᵀ * s.Zu N ω)ᵀ) := by
  rw [s.hE N ω, Matrix.transpose_smul, Matrix.smul_mul, Matrix.transpose_mul,
    Matrix.transpose_transpose]

/-- The `d × r` matrix `Q = V Q_C Λ^{1/2} + Eᵀ U`. Column `j` is `√λ_j (V q_j) + g_j`. -/
noncomputable def qmatR (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (ns N)) (Fin r) ℝ) : Matrix (Fin (d N)) (Fin r) ℝ :=
  s.spikeMat N * Matrix.diagonal (fun j => Real.sqrt (s.coreEig j)) + (s.E N ω)ᵀ * U

theorem qmatR_apply (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (ns N)) (Fin r) ℝ) (k : Fin (d N)) (j : Fin r) :
    s.qmatR N ω U k j
      = Real.sqrt (s.coreEig j) * WithLp.ofLp (s.spikeVec j N) k
        + ((s.E N ω)ᵀ * U) k j := by
  rw [qmatR, Matrix.add_apply, Matrix.mul_diagonal, s.spikeMat_apply N k j, mul_comm]

/-- The frame complement of the noise, `E⊥ = E - U Uᵀ E`. -/
noncomputable def stackEperp (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (ns N)) (Fin r) ℝ) : Matrix (Fin (ns N)) (Fin (d N)) ℝ :=
  s.E N ω - U * (Uᵀ * s.E N ω)

/-- The rank-`r` Wishart block `W₀ = E⊥ᵀ E⊥` of the split. -/
noncomputable def rankRW0 (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (ns N)) (Fin r) ℝ) : Matrix (Fin (d N)) (Fin (d N)) ℝ :=
  (s.stackEperp N ω U)ᵀ * s.stackEperp N ω U

theorem isHermitian_rankRW0 (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (ns N)) (Fin r) ℝ) : (s.rankRW0 N ω U).IsHermitian :=
  isHermitian_transpose_mul_self _

/-- The block `W₀` of the split is positive semidefinite. -/
theorem eigenvalues_rankRW0_nonneg (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    (U : Matrix (Fin (ns N)) (Fin r) ℝ) (a : Fin (d N)) :
    0 ≤ (s.isHermitian_rankRW0 N ω U).eigenvalues a := by
  have hpsd : (s.rankRW0 N ω U).PosSemidef := by
    change ((s.stackEperp N ω U)ᵀ * s.stackEperp N ω U).PosSemidef
    simpa using Matrix.posSemidef_conjTranspose_mul_self (s.stackEperp N ω U)
  exact hpsd.eigenvalues_nonneg a

/-- The table as frame times coefficients plus complement, for any `U` with the factor
property. -/
theorem X_eq_frame_add (s : RankRStack μ ns d r) (N : ℕ) (ω : Ω N)
    {U : Matrix (Fin (ns N)) (Fin r) ℝ}
    (hsig : s.signalPart N
      = U * (s.spikeMat N * Matrix.diagonal fun j => Real.sqrt (s.coreEig j))ᵀ) :
    s.X N ω = U * (s.qmatR N ω U)ᵀ + s.stackEperp N ω U := by
  have hT : ((s.E N ω)ᵀ * U)ᵀ = Uᵀ * s.E N ω := by
    rw [Matrix.transpose_mul, Matrix.transpose_transpose]
  rw [s.X_eq N ω, hsig, qmatR, Matrix.transpose_add, Matrix.mul_add, hT, stackEperp]
  abel

/-! ### 5. Measurability -/

theorem measurable_E (s : RankRStack μ ns d r) (N : ℕ) : Measurable (s.E N) := by
  refine measurable_pi_lambda _ fun q => measurable_pi_lambda _ fun l => ?_
  have h : (fun ω => s.E N ω q l)
      = fun ω => (Real.sqrt (d N))⁻¹ * s.Zu N ω q l := by
    funext ω
    rw [s.hE N ω]
    rfl
  rw [h]
  exact ((measurable_pi_apply l).comp
    ((measurable_pi_apply q).comp (s.hZmeas N))).const_mul _

theorem measurable_stackEperp (s : RankRStack μ ns d r) (N : ℕ)
    (U : Matrix (Fin (ns N)) (Fin r) ℝ) :
    Measurable (fun ω => s.stackEperp N ω U) := by
  have hE : ∀ (q : Fin (ns N)) (l : Fin (d N)), Measurable fun ω => s.E N ω q l :=
    fun q l => ((measurable_pi_apply l).comp ((measurable_pi_apply q).comp (s.measurable_E N)))
  refine measurable_pi_lambda _ fun q => measurable_pi_lambda _ fun l => ?_
  have h : (fun ω => s.stackEperp N ω U q l)
      = fun ω => s.E N ω q l
        - ∑ a : Fin r, U q a * ∑ b, U b a * s.E N ω b l := by
    funext ω
    simp only [stackEperp, Matrix.sub_apply, Matrix.mul_apply, Matrix.transpose_apply]
  rw [h]
  refine (hE q l).sub ?_
  refine Finset.measurable_sum _ fun a _ => ?_
  exact (Finset.measurable_sum _ fun b _ => (hE b l).const_mul (U b a)).const_mul (U q a)

theorem measurableSet_lamMax_rankRW0_le (s : RankRStack μ ns d r) (N : ℕ)
    (U : Matrix (Fin (ns N)) (Fin r) ℝ) (x : ℝ) :
    MeasurableSet {ω | lamMax (s.rankRW0 N ω U) (s.isHermitian_rankRW0 N ω U) ≤ x} := by
  have h : {ω | lamMax (s.rankRW0 N ω U) (s.isHermitian_rankRW0 N ω U) ≤ x}
      = (fun ω => gramLamMax (s.stackEperp N ω U)) ⁻¹' Set.Iic x := rfl
  rw [h]
  exact (measurable_gramLamMax.comp (s.measurable_stackEperp N U)) measurableSet_Iic

end RankRStack

end StackedSVD
