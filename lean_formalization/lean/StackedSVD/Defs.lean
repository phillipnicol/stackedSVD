/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import Mathlib

/-!
# Definitions v2: rank-one spiked model, spectral projectors, overlap

STATUS 2026-08-29: compiled on Mathlib v4.33.0 (lake env lean, server).
See `notes/archive/prop_single_table.md` (modeling choices) and `notes/archive/SINGLE_TABLE_PLAN.md`
§3.

Changes from v1 (all from the audit, user OK 2026-08-29):
* one probability space `(Ω N, μ N)` per `N`; `TendstoInProb` is a tail-probability statement;
* rows are `Fin (n N)`, vectors are `EuclideanSpace`, dimensions are positive;
* noise is stated as a law (`HasLaw` against the canonical Gaussian matrix measure), so that
  Layer 2 proves everything about the canonical measure and transfers for free;
* `specProj A S` is the basis-free spectral projector (eigenvalues in `S`); `topProj` is the
  case `S = {lamMax A}`; `overlap X w = ‖topProj (Xᵀ X) w‖²`.

Added on the server (2026-08-29) to make the file elaborate: `instMeasurableSpaceMatrix`,
the product sigma-algebra on `Matrix m n α`. Mathlib v4.33.0 has no such instance.

`isHermitian_mul_transpose_self` is added here (F28, 2026-09-08), the dedup of the
copies formerly in `SVDStack/Defs.lean` and `RMT/Simplicity.lean`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped InnerProductSpace Matrix

namespace StackedSVD

/-! ### Convergence in probability, one space per `N` -/

/-- Convergence in probability to a constant along `N → ∞`, with a probability space
`(Ω N, μ N)` for each `N`. -/
def TendstoInProb {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] (μ : ∀ N, Measure (Ω N))
    (f : ∀ N, Ω N → ℝ) (a : ℝ) : Prop :=
  ∀ ε > 0, Tendsto (fun N => μ N {ω | ε ≤ |f N ω - a|}) atTop (𝓝 0)

/-! ### Scalars -/

/-- `β²` of `prop:single_table`. The paper's threshold `θ ≥ c^{1/4}` equals `θ⁴ ≥ c` for
`θ ≥ 0`; at equality both branches give `0`, so `>` is used. -/
noncomputable def betaSq (θ c : ℝ) : ℝ :=
  if θ ^ 4 > c then (θ ^ 4 - c) / (θ ^ 4 + θ ^ 2) else 0

/-- Bulk edge of `W = Eᵀ E` under the `1/d` variance scaling. -/
noncomputable def bulkEdge (c : ℝ) : ℝ := (1 + Real.sqrt c) ^ 2

/-- Limit of the top eigenvalue of `Xᵀ X`: the outlier `θ² + 1 + c + c/θ²` above threshold,
the bulk edge `bulkEdge c = (1 + √c)²` otherwise. (Both equal `(1+√c)²` at `θ⁴ = c`.) -/
noncomputable def rhoSq (θ c : ℝ) : ℝ :=
  if θ ^ 4 > c then θ ^ 2 + 1 + c + c / θ ^ 2 else bulkEdge c

/-! ### Spectral projectors and overlap (basis-free) -/

section Spectral

variable {d : ℕ}

/-- The symmetric operator of a real symmetric matrix on `EuclideanSpace ℝ (Fin d)`. -/
noncomputable abbrev toOp (A : Matrix (Fin d) (Fin d) ℝ) :
    EuclideanSpace ℝ (Fin d) →ₗ[ℝ] EuclideanSpace ℝ (Fin d) :=
  Matrix.toEuclideanLin A

/-- Largest eigenvalue of a real symmetric matrix (`eigenvalues₀` is antitone, index `0`).
Junk value `0` when `d = 0`. -/
noncomputable def lamMax (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian) : ℝ :=
  if h : 0 < d then hA.eigenvalues₀ ⟨0, by simpa using h⟩ else 0

/-- Sum of the eigenspaces of `A` with eigenvalue in `S`. -/
noncomputable def specSpace (A : Matrix (Fin d) (Fin d) ℝ) (S : Set ℝ) :
    Submodule ℝ (EuclideanSpace ℝ (Fin d)) :=
  ⨆ (t : ℝ) (_ : t ∈ S), Module.End.eigenspace (toOp A) t

/-- Orthogonal projector (as an operator `E → E`) onto `specSpace A S`. Basis-free.
API name to confirm on the pinned Mathlib: `Submodule.starProjection` (v4.33) or
`(orthogonalProjection K).subtype ∘L ...`; `HasOrthogonalProjection` holds for
finite-dimensional submodules. -/
noncomputable def specProj (A : Matrix (Fin d) (Fin d) ℝ) (S : Set ℝ) :
    EuclideanSpace ℝ (Fin d) →L[ℝ] EuclideanSpace ℝ (Fin d) :=
  (specSpace A S).starProjection

/-- Top eigenspace and its projector. -/
noncomputable def topSpace (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian) :
    Submodule ℝ (EuclideanSpace ℝ (Fin d)) :=
  specSpace A {lamMax A hA}

noncomputable def topProj (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian) :
    EuclideanSpace ℝ (Fin d) →L[ℝ] EuclideanSpace ℝ (Fin d) :=
  specProj A {lamMax A hA}

/-- The top eigenvalue is simple. -/
def TopSimple (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian) : Prop :=
  Module.finrank ℝ (topSpace A hA) = 1

/-- The Gram matrix `Xᵀ X` is Hermitian (real symmetric). -/
theorem isHermitian_transpose_mul_self {n : ℕ} (X : Matrix (Fin n) (Fin d) ℝ) :
    (Xᵀ * X).IsHermitian := by
  simpa using Matrix.isHermitian_conjTranspose_mul_self X

/-- Companion of `isHermitian_transpose_mul_self` of `Defs.lean`, for `A Aᵀ`. -/
theorem isHermitian_mul_transpose_self {M d : ℕ} (A : Matrix (Fin M) (Fin d) ℝ) :
    (A * Aᵀ).IsHermitian := by
  simpa using Matrix.isHermitian_mul_conjTranspose_self A

/-- Squared overlap of the top right singular subspace of `X` with `w`: `‖P w‖²` with
`P = topProj (Xᵀ X)`. Equals `⟨v̂, w⟩²` for a unit top eigenvector `v̂` when the top
eigenvalue is simple (`overlap_eq_inner_sq`, to prove). -/
noncomputable def overlap {n : ℕ} (X : Matrix (Fin n) (Fin d) ℝ)
    (w : EuclideanSpace ℝ (Fin d)) : ℝ :=
  ‖topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) w‖ ^ 2

/-- Top eigenvalue of the Gram matrix. -/
noncomputable def gramLamMax {n : ℕ} (X : Matrix (Fin n) (Fin d) ℝ) : ℝ :=
  lamMax (Xᵀ * X) (isHermitian_transpose_mul_self X)

end Spectral

/-! ### Canonical Gaussian matrix measure -/

/-- Product sigma-algebra on matrices. Mathlib v4.33.0 gives `Matrix m n α` no
`MeasurableSpace` instance: `Matrix` is a `def`, so every instance is declared by hand
(compare `Matrix.decidableEq`, `Fintype (Matrix m n α)` in `Mathlib/Data/Matrix/Basic.lean`).
`inferInstanceAs` keeps this definitionally equal to `MeasurableSpace.pi`, which is the
sigma-algebra that `Measure.pi` builds. -/
instance instMeasurableSpaceMatrix {m n α : Type*} [MeasurableSpace α] :
    MeasurableSpace (Matrix m n α) :=
  inferInstanceAs (MeasurableSpace (m → n → α))

/-- Law of an `n × d` matrix with i.i.d. `N(0,1)` entries (product of `gaussianReal 0 1`). -/
noncomputable def gaussianMatrix (n d : ℕ) : Measure (Matrix (Fin n) (Fin d) ℝ) :=
  Measure.pi fun _ : Fin n => Measure.pi fun _ : Fin d => gaussianReal 0 1

/-! ### The spiked model -/

/-- Rank-one spiked model on one probability space per `N`: `X N = θ u_N v_Nᵀ + d_N^{-1/2} Z N`.
`Z` is the unscaled noise (entries `N(0,1)` under `GaussianNoise`). -/
structure SpikedModel {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    (μ : ∀ N, Measure (Ω N)) (n d : ℕ → ℕ) where
  /-- signal strength -/
  θ : ℝ
  /-- left singular vector, deterministic -/
  u : (N : ℕ) → EuclideanSpace ℝ (Fin (n N))
  /-- right singular vector, deterministic, shared across tables -/
  v : (N : ℕ) → EuclideanSpace ℝ (Fin (d N))
  /-- unscaled noise -/
  Z : (N : ℕ) → Ω N → Matrix (Fin (n N)) (Fin (d N)) ℝ
  hθ : 0 ≤ θ
  hn : ∀ N, 0 < n N
  hd : ∀ N, 0 < d N
  hu : ∀ N, ‖u N‖ = 1
  hv : ∀ N, ‖v N‖ = 1
  hZ : ∀ N, Measurable (Z N)

namespace SpikedModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)} {n d : ℕ → ℕ}

/-- Scaled noise `E_N = d_N^{-1/2} Z_N`. -/
noncomputable def E (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (n N)) (Fin (d N)) ℝ :=
  (Real.sqrt (d N))⁻¹ • m.Z N ω

/-- `X_N ω = θ u_N v_Nᵀ + E_N ω`. -/
noncomputable def X (m : SpikedModel μ n d) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (n N)) (Fin (d N)) ℝ :=
  m.θ • Matrix.vecMulVec (WithLp.ofLp (m.u N)) (WithLp.ofLp (m.v N)) + m.E N ω

/-- `assum:general_noise`, Gaussian case, as a law: `Z_N ~ gaussianMatrix`. -/
def GaussianNoise (m : SpikedModel μ n d) : Prop :=
  ∀ N, HasLaw (m.Z N) (gaussianMatrix (n N) (d N)) (μ N)

/-- `eq:RMT_limit` for one table: `n → ∞`, `d → ∞`, `n / d → c`. -/
def Regime (_m : SpikedModel μ n d) (c : ℝ) : Prop :=
  Tendsto n atTop atTop ∧ Tendsto d atTop atTop ∧
    Tendsto (fun N => (n N : ℝ) / d N) atTop (𝓝 c)

end SpikedModel

/-! ### Several tables on the same spaces -/

/-- `M` tables with a shared `v`, all on `(Ω N, μ N)`. Cross-table independence is a
separate predicate `JointGaussianNoise`. -/
structure MultiTableModel {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    (μ : ∀ N, Measure (Ω N)) (M : ℕ) (n : Fin M → ℕ → ℕ) (d : ℕ → ℕ) where
  tbl : (i : Fin M) → SpikedModel μ (n i) d
  hv : ∀ i j N, (tbl i).v N = (tbl j).v N

namespace MultiTableModel

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ}

/-- Joint law of all noise matrices at level `N`: independent Gaussian tables.
This is `assum:general_noise` (Gaussian) for all tables at once, including independence
across tables, which the single-table predicate cannot express. -/
def JointGaussianNoise (m : MultiTableModel μ M n d) : Prop :=
  ∀ N, HasLaw (fun ω i => (m.tbl i).Z N ω)
    (Measure.pi fun i : Fin M => gaussianMatrix (n i N) (d N)) (μ N)

end MultiTableModel

end StackedSVD
