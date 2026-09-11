/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.Split
import StackedSVD.RankR.RMT.Split
import StackedSVD.RankR.StackGamma
import StackedSVD.RankR.SubspaceGStack

/-!
# Stage E1: the column split of the weighted rank-`r` stack

Stage E1 of `notes/archive/rankr_TrackE_plan.md` (section 2.1, steps 1 and 2 of section 2.3, and the
E1 row of section 3). The rank-one mirror is `RMT/Het/Split.lean`, which runs the same split
on a `MultiTableModel` with one spike. Here the `r` columns of `V` replace the single vector
`v`, and a product `Q Qᵀ` of an `n_tot × r` matrix replaces `vecMulVec q q`.

The weighted stack of an `UnalignedModelR` is

```
X_w = [w_1 X_1; ...; w_M X_M] = Σ^{1/2} (A Vᵀ + E_stack),   Σ^{1/2} = diag(w_i I_{n_i}),
```

with `A = signalFactorG` the `n_tot × r` signal factor (block row `i` is `U_i Θ_i R_iᵀ`) and
`E_stack` the **unweighted** stacked noise. This file splits the **columns** off `V`:

```
E⊥ = E_stack (1 - V Vᵀ),   Ũ = Σ^{1/2} A,   Q = Ũ + Σ^{1/2} E_stack V,
W₀' = Σ^{1/2} E⊥ E⊥ᵀ Σ^{1/2},
```

and proves `X_w X_wᵀ = W₀' + Q Qᵀ` with `E⊥ V = 0`. The update is an honest rank-`r` positive
update of a block that is independent of `Q`, which is what the rank-`r` deterministic layer
of `RankR/RMT/OutliersG.lean` consumes.

The row split of `RankR/RMT/Stack.lean` is not available here, for the reason its own header
gives at rank one: `(1 - U Uᵀ) Σ^{1/2} E` is independent of `Uᵀ Σ^{1/2} E` only when the
column span of `U` is invariant under `Σ`, which the weighted stack does not give.

Two modeling choices are inherited from the rank-one file:

1. `SigmaHalfR` carries the **signed** weight `w_i`, not `|w_i|`. It is then exactly the
   matrix that turns the unweighted stack into the weighted one
   (`stackXW_eq_SigmaHalf_mul`). The paper's `Σ^{1/2} = diag(|w_i|)` differs by a per-block
   sign that the sign symmetry of the noise law absorbs.
2. `EperpHetR` keeps all `d N` columns and kills the column span of `V`
   (`EperpHetR_mul_V`), rather than living on a `d - r` dimensional complement. The `d - r`
   block appears only inside `exists_block_hasLaw_hetR`, where the law is stated.

Paper: `main_paper.tex:2294` (`eq:rank_r_model`), `:2325` (`eq:stacksvd_appXstack`).

No `sorry`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-! ### 1. The five split objects -/

/-- `Σ^{1/2} = diag(w_{blk(q)})`: the weight of the table that owns row `q`, through the same
`finSigmaFinEquiv` row order as `stackXG`. Signed, see choice 1 of the header. The rank-one
mirror is `MultiTableModel.SigmaHalf` (`RMT/Het/Split.lean:168`). -/
noncomputable def SigmaHalfR (_m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ) :
    Matrix (Fin (∑ i, n i N)) (Fin (∑ i, n i N)) ℝ :=
  Matrix.diagonal fun q => w (finSigmaFinEquiv.symm q).1

/-- `Ũ = Σ^{1/2} A`, the weighted signal factor (`n_tot × r`). At `r = 1` this is
`MultiTableModel.u0Het` (`RMT/Het/Split.lean:174`) read as a one-column matrix. -/
noncomputable def UtildeR (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ) :
    Matrix (Fin (∑ i, n i N)) (Fin r) ℝ :=
  m.SigmaHalfR w N * m.signalFactorG N

/-- `E⊥ = E_stack (1 - V Vᵀ)`: the unweighted stacked noise with its `V` component removed.
The rank-one mirror is `MultiTableModel.EperpHet` (`RMT/Het/Split.lean:184`). -/
noncomputable def EperpHetR (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    Matrix (Fin (∑ i, n i N)) (Fin (d N)) ℝ :=
  m.stackEG N ω * (1 - m.V N * (m.V N)ᵀ)

/-- `Q = Ũ + Σ^{1/2} E_stack V` (`n_tot × r`): column `k` is `X_w v_k`
(`stackXW_mulVec_colVecG`). The rank-one mirror is `MultiTableModel.qHet`
(`RMT/Het/Split.lean:189`). -/
noncomputable def QmatHetR (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ :=
  m.UtildeR w N + m.SigmaHalfR w N * (m.stackEG N ω * m.V N)

/-- `W₀' = Σ^{1/2} E⊥ E⊥ᵀ Σ^{1/2}`, the noise Gram matrix on the `n` side. The rank-one
mirror is `MultiTableModel.W0het` (`RMT/Het/Split.lean:194`). -/
noncomputable def W0hetR (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : Matrix (Fin (∑ i, n i N)) (Fin (∑ i, n i N)) ℝ :=
  m.SigmaHalfR w N * m.EperpHetR N ω * (m.EperpHetR N ω)ᵀ * m.SigmaHalfR w N

/-! ### 2. Deterministic algebra -/

theorem transpose_SigmaHalfR (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ) :
    (m.SigmaHalfR w N)ᵀ = m.SigmaHalfR w N :=
  Matrix.diagonal_transpose _

/-- The rank-`r` `Σ^{1/2}` is the rank-one one of the shell `toMultiTableShell`. Both are the
same diagonal matrix, so the identity is `rfl`; it is here so that a later stage may quote a
rank-one lemma of `RMT/Het/Split.lean` on this object. -/
theorem SigmaHalfR_eq_shell [NeZero M] (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ)
    (N : ℕ) : m.SigmaHalfR w N = m.toMultiTableShell.SigmaHalf w N := rfl

/-- `X_w = Σ^{1/2} X_stack`: the weighted stack is the unweighted one with the rows of block
`i` scaled by `w_i`. The rank-one mirror is `stackZW_eq_SigmaHalf_mul`
(`RMT/Het/Split.lean:215`). -/
theorem stackXW_eq_SigmaHalf_mul (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : m.stackXW w N ω = m.SigmaHalfR w N * m.stackXG N ω := by
  ext q k
  rw [SigmaHalfR, Matrix.diagonal_mul, stackXW_apply', stackXG_apply']

/-- `E⊥ V = 0`: the column analogue of `EperpHet_mulVec_v` (`RMT/Het/Split.lean:235`) at `r`
columns. The only model field it reads is `m.hV N`. -/
theorem EperpHetR_mul_V (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    m.EperpHetR N ω * m.V N = 0 := by
  rw [EperpHetR, Matrix.mul_assoc, Matrix.sub_mul, Matrix.one_mul, Matrix.mul_assoc, m.hV N,
    Matrix.mul_one, sub_self, Matrix.mul_zero]

/-- `X_w = Q Vᵀ + Σ^{1/2} E⊥`, the rank-`r` twin of `stackW_X_eq_vecMulVec_add`. -/
theorem stackXW_eq_add (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) :
    m.stackXW w N ω
      = m.QmatHetR w N ω * (m.V N)ᵀ + m.SigmaHalfR w N * m.EperpHetR N ω := by
  rw [m.stackXW_eq_SigmaHalf_mul w N ω, m.stackX_eqG N ω, signalPartG, QmatHetR, UtildeR,
    EperpHetR]
  simp only [Matrix.mul_add, Matrix.add_mul, Matrix.mul_sub, Matrix.mul_one, Matrix.mul_assoc]
  abel

/-- `X_w V = Q`: the `r` columns of `Q` are the images of the columns of `V`. -/
theorem stackXW_mul_V (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    m.stackXW w N ω * m.V N = m.QmatHetR w N ω := by
  rw [m.stackXW_eq_add w N ω, Matrix.add_mul, Matrix.mul_assoc, m.hV N, Matrix.mul_one,
    Matrix.mul_assoc, m.EperpHetR_mul_V N ω, Matrix.mul_zero, add_zero]

/-- `X_w v_k = Q e_k`: the vector the rank-`r` duality of stage E2 reads (step 6 of section
2.3 of the plan). -/
theorem stackXW_mulVec_colVecG (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) (k : Fin r) :
    m.stackXW w N ω *ᵥ WithLp.ofLp (m.colVecG N k) = fun q => m.QmatHetR w N ω q k := by
  have hv : WithLp.ofLp (m.colVecG N k) = fun l => m.V N l k := rfl
  rw [hv]
  funext q
  rw [← m.stackXW_mul_V w N ω, Matrix.mul_apply]
  rfl

/-- `W₀' = Y Yᵀ` with `Y = Σ^{1/2} E⊥`. The rank-one mirror is `W0het_eq_mul_transpose`
(`RMT/Het/R5het.lean:376`). -/
theorem W0hetR_eq_mul_transpose (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) :
    m.W0hetR w N ω
      = (m.SigmaHalfR w N * m.EperpHetR N ω) * (m.SigmaHalfR w N * m.EperpHetR N ω)ᵀ := by
  rw [W0hetR, Matrix.transpose_mul, m.transpose_SigmaHalfR w N, Matrix.mul_assoc,
    Matrix.mul_assoc]

/-- **Step 2 of the plan.** `X_w X_wᵀ = W₀' + Q Qᵀ`. Pure algebra: it reads `Vᵀ V = 1` and
`E⊥ V = 0` only. The rank-one mirror is `gram_eq_het` (`RMT/Het/Split.lean:282`). -/
theorem gram_eq_hetR (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ) (ω : Ω N) :
    m.stackXW w N ω * (m.stackXW w N ω)ᵀ
      = m.W0hetR w N ω + m.QmatHetR w N ω * (m.QmatHetR w N ω)ᵀ := by
  have hGV : m.SigmaHalfR w N * m.EperpHetR N ω * m.V N = 0 := by
    rw [Matrix.mul_assoc, m.EperpHetR_mul_V N ω, Matrix.mul_zero]
  rw [m.stackXW_eq_add w N ω, m.W0hetR_eq_mul_transpose w N ω]
  set Q := m.QmatHetR w N ω with hQ
  set G := m.SigmaHalfR w N * m.EperpHetR N ω with hG
  rw [Matrix.transpose_add, Matrix.transpose_mul, Matrix.transpose_transpose, Matrix.add_mul,
    Matrix.mul_add, Matrix.mul_add]
  have e1 : Q * (m.V N)ᵀ * (m.V N * Qᵀ) = Q * Qᵀ := by
    rw [← Matrix.mul_assoc, Matrix.mul_assoc Q, m.hV N, Matrix.mul_one]
  have e2 : Q * (m.V N)ᵀ * Gᵀ = 0 := by
    rw [Matrix.mul_assoc, ← Matrix.transpose_mul, hGV, Matrix.transpose_zero,
      Matrix.mul_zero]
  have e3 : G * (m.V N * Qᵀ) = 0 := by
    rw [← Matrix.mul_assoc, hGV, Matrix.zero_mul]
  rw [e1, e2, e3]
  abel

theorem isHermitian_W0hetR (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : (m.W0hetR w N ω).IsHermitian := by
  rw [m.W0hetR_eq_mul_transpose w N ω]
  simpa using
    Matrix.isHermitian_mul_conjTranspose_self (m.SigmaHalfR w N * m.EperpHetR N ω)

theorem posSemidef_W0hetR (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : (m.W0hetR w N ω).PosSemidef := by
  rw [m.W0hetR_eq_mul_transpose w N ω]
  simpa using
    Matrix.posSemidef_self_mul_conjTranspose (m.SigmaHalfR w N * m.EperpHetR N ω)

/-! ### 3. The block law -/

/-- `E⊥ E⊥ᵀ = d⁻¹ (Z Zᵀ - (Z V)(Z V)ᵀ)` on the unscaled stacked noise. The rank-one mirror is
`EperpHet_gram_eq_smul` (`RMT/Het/Split.lean:344`). No positivity of `d N` is needed: at
`d N = 0` both sides are `0` and `(√0)⁻¹ (√0)⁻¹ = 0 = (0 : ℝ)⁻¹`. -/
theorem EperpHetR_gram_eq_smul (m : UnalignedModelR μ M n d r rk) (N : ℕ) (ω : Ω N) :
    m.EperpHetR N ω * (m.EperpHetR N ω)ᵀ
      = ((d N : ℝ))⁻¹ • (m.stackZG N ω * (m.stackZG N ω)ᵀ
          - (m.stackZG N ω * m.V N) * (m.stackZG N ω * m.V N)ᵀ) := by
  have hsq : (Real.sqrt (d N))⁻¹ * (Real.sqrt (d N))⁻¹ = ((d N : ℝ))⁻¹ := by
    rw [← mul_inv, Real.mul_self_sqrt (Nat.cast_nonneg _)]
  have hVV : (m.V N * (m.V N)ᵀ) * (m.V N * (m.V N)ᵀ) = m.V N * (m.V N)ᵀ := by
    rw [Matrix.mul_assoc, ← Matrix.mul_assoc (m.V N)ᵀ, m.hV N, Matrix.one_mul]
  have hPt : ((1 : Matrix (Fin (d N)) (Fin (d N)) ℝ) - m.V N * (m.V N)ᵀ)ᵀ
      = 1 - m.V N * (m.V N)ᵀ := by
    rw [Matrix.transpose_sub, Matrix.transpose_one, Matrix.transpose_mul,
      Matrix.transpose_transpose]
  have hPP : ((1 : Matrix (Fin (d N)) (Fin (d N)) ℝ) - m.V N * (m.V N)ᵀ)
      * (1 - m.V N * (m.V N)ᵀ) = 1 - m.V N * (m.V N)ᵀ := by
    rw [Matrix.sub_mul, Matrix.one_mul, Matrix.mul_sub, Matrix.mul_one, hVV]
    abel
  have hE : m.EperpHetR N ω
      = (Real.sqrt (d N))⁻¹ • (m.stackZG N ω * (1 - m.V N * (m.V N)ᵀ)) := by
    rw [EperpHetR, m.stackE_eqG N ω, Matrix.smul_mul]
  rw [hE, Matrix.transpose_smul, Matrix.smul_mul, Matrix.mul_smul, smul_smul, hsq]
  congr 1
  rw [Matrix.transpose_mul, hPt, ← Matrix.mul_assoc, Matrix.mul_assoc (m.stackZG N ω), hPP,
    Matrix.mul_sub, Matrix.mul_one, Matrix.sub_mul, Matrix.transpose_mul]
  simp only [Matrix.mul_assoc]

/-- `W₀'` written on a block `B` with `E⊥ E⊥ᵀ = d⁻¹ B Bᵀ`. Deterministic bookkeeping: it turns
the block identity of `exists_block_hasLaw_hetR` into the form the rank-`r` deterministic
layer consumes. The rank-one mirror is `W0het_eq_of_block` (`RMT/Het/Split.lean:357`). -/
theorem W0hetR_eq_of_block (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) {p : ℕ} (B : Matrix (Fin (∑ i, n i N)) (Fin p) ℝ)
    (hB : m.EperpHetR N ω * (m.EperpHetR N ω)ᵀ = ((d N : ℝ))⁻¹ • (B * Bᵀ)) :
    m.W0hetR w N ω
      = ((d N : ℝ))⁻¹ • ((m.SigmaHalfR w N * B) * (m.SigmaHalfR w N * B)ᵀ) := by
  have h1 : m.W0hetR w N ω
      = m.SigmaHalfR w N * (m.EperpHetR N ω * (m.EperpHetR N ω)ᵀ) * m.SigmaHalfR w N := by
    simp only [W0hetR, Matrix.mul_assoc]
  rw [h1, hB, Matrix.mul_smul, Matrix.smul_mul]
  congr 1
  rw [Matrix.transpose_mul, m.transpose_SigmaHalfR w N]
  simp only [Matrix.mul_assoc]

/-- **The block law, step 3 of the plan.** For `d N = p + r` there is an `n_tot × p` block `B`
with `E⊥ E⊥ᵀ = d⁻¹ B Bᵀ`, and the pair `(Z_stack V, B)` has a product Gaussian law, which
carries both marginals and their independence at once. The rank-one mirror is
`exists_block_hasLaw_het` (`RMT/Het/Split.lean:390`).

Route (section 2.3 step 1 of the plan). The column split is the row split of the transpose,
so `exists_frameBlock'` (`RankR/RMT/Split.lean:446`) runs at `nn = d N`, `rr = r`,
`d := ∑ i, n i N`, `U = V N`, on `Z := Z_stackᵀ`:

1. `hasLaw_stackZG` gives `Z_stack ~ gaussianMatrix (∑ n_i) (d N)`.
2. `measurePreserving_transpose` moves it to `Z_stackᵀ ~ gaussianMatrix (d N) (∑ n_i)`.
3. `exists_frameBlock'` completes the columns of `V` to an orthonormal basis of `ℝ^{d N}`,
   rotates, and splits the first `r` rows off the block of the last `p` rows. Its first
   component is `Zᵀ U = Z_stack V` already. Put `B := (Yfun Z_stackᵀ)ᵀ`.
4. A final `Prod.map id transpose` turns `gaussianMatrix p (∑ n_i)` into
   `gaussianMatrix (∑ n_i) p`.
5. The Gram half of step 3 transposes to `B Bᵀ = Z Zᵀ - (Z V)(Z V)ᵀ`, and
   `EperpHetR_gram_eq_smul` closes the first conjunct.
6. `p = 0` needs no special case: `B` is the empty matrix and both sides of the first
   conjunct are `0`. -/
theorem exists_block_hasLaw_hetR [NeZero M] (m : UnalignedModelR μ M n d r rk)
    (hG : m.JointGaussianNoise) (N : ℕ) {p : ℕ} (hp : d N = p + r) :
    ∃ B : Ω N → Matrix (Fin (∑ i, n i N)) (Fin p) ℝ,
      (∀ ω, m.EperpHetR N ω * (m.EperpHetR N ω)ᵀ = ((d N : ℝ))⁻¹ • (B ω * (B ω)ᵀ)) ∧
        HasLaw (fun ω => (m.stackZG N ω * m.V N, B ω))
          ((gaussianMatrix (∑ i, n i N) r).prod (gaussianMatrix (∑ i, n i N) p)) (μ N) := by
  obtain ⟨Yfun, hgram, hmp⟩ :=
    exists_frameBlock' (d := ∑ i, n i N) (hp.trans (add_comm p r)) (m.V N) (m.hV N)
  refine ⟨fun ω => (Yfun ((m.stackZG N ω)ᵀ))ᵀ, ?_, ?_⟩
  · intro ω
    have h := hgram ((m.stackZG N ω)ᵀ)
    rw [Matrix.transpose_transpose] at h
    rw [m.EperpHetR_gram_eq_smul N ω, Matrix.transpose_transpose, h]
  · have hZt : HasLaw (fun ω => (m.stackZG N ω)ᵀ)
        (gaussianMatrix (d N) (∑ i, n i N)) (μ N) :=
      (measurePreserving_transpose (∑ i, n i N) (d N)).fun_comp_hasLaw (m.hasLaw_stackZG hG N)
    have hprod : MeasurePreserving
        (Prod.map (id : Matrix (Fin (∑ i, n i N)) (Fin r) ℝ →
            Matrix (Fin (∑ i, n i N)) (Fin r) ℝ)
          (Matrix.transpose :
            Matrix (Fin p) (Fin (∑ i, n i N)) ℝ → Matrix (Fin (∑ i, n i N)) (Fin p) ℝ))
        ((gaussianMatrix (∑ i, n i N) r).prod (gaussianMatrix p (∑ i, n i N)))
        ((gaussianMatrix (∑ i, n i N) r).prod (gaussianMatrix (∑ i, n i N) p)) :=
      (MeasurePreserving.id _).prod (measurePreserving_transpose p (∑ i, n i N))
    have h := (hprod.comp hmp).fun_comp_hasLaw hZt
    simpa [Function.comp_def] using h

/-! ### 4. Measurability

The facts stages E3, E4, E5 and E7 read. `stackXW` is affine in the unscaled noise entry by
entry, `QmatHetR` is `stackXW * V`, and `W₀'` is a product of two blocks that are linear in
the noise. The edge event of stage E3 is `measurableSet_lamMax_W0hetR_le`. -/

/-- The weighted stack is measurable, entry by entry from `SpikedModelR.hZ`. -/
theorem measurable_stackXW (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ) :
    Measurable (m.stackXW w N) := by
  refine measurable_pi_lambda _ fun q => measurable_pi_lambda _ fun k => ?_
  simp only [stackXW_apply', SpikedModelR.X, SpikedModelR.E, Matrix.add_apply,
    Matrix.smul_apply, smul_eq_mul]
  exact measurable_const.mul (measurable_const.add (measurable_const.mul
    ((measurable_pi_apply k).comp ((measurable_pi_apply _).comp
      ((m.tbl (finSigmaFinEquiv.symm q).1).hZ N)))))

/-- `Q` is measurable, through `stackXW_mul_V`. -/
theorem measurable_QmatHetR (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ) :
    Measurable (m.QmatHetR w N) := by
  have h : m.QmatHetR w N = fun ω => m.stackXW w N ω * m.V N := by
    funext ω
    rw [m.stackXW_mul_V w N ω]
  rw [h]
  refine measurable_pi_lambda _ fun q => measurable_pi_lambda _ fun k => ?_
  simp only [Matrix.mul_apply]
  refine Finset.measurable_sum _ fun l _ => ?_
  exact (((measurable_pi_apply l).comp
    ((measurable_pi_apply q).comp (m.measurable_stackXW w N))).mul_const _)


/-- The stacked scaled noise is measurable, from `measurable_stackZG`
(`RankR/SubspaceGStack.lean:89`). -/
theorem measurable_stackEG (m : UnalignedModelR μ M n d r rk) (N : ℕ) :
    Measurable (m.stackEG N) := by
  refine measurable_pi_lambda _ fun q => measurable_pi_lambda _ fun l => ?_
  have h : (fun ω => m.stackEG N ω q l)
      = fun ω => (Real.sqrt (d N))⁻¹ * m.stackZG N ω q l := by
    funext ω
    rw [m.stackE_eqG N ω, Matrix.smul_apply, smul_eq_mul]
  rw [h]
  exact ((measurable_pi_apply l).comp
    ((measurable_pi_apply q).comp (m.measurable_stackZG N))).const_mul _

/-- `E⊥` is measurable, entry by entry. The rank-one mirror is `measurable_EperpHet`
(`RMT/Het/R5het.lean:382`). -/
theorem measurable_EperpHetR (m : UnalignedModelR μ M n d r rk) (N : ℕ) :
    Measurable (m.EperpHetR N) := by
  refine measurable_pi_lambda _ fun q => measurable_pi_lambda _ fun j => ?_
  have h : (fun ω => m.EperpHetR N ω q j)
      = fun ω => ∑ l, m.stackEG N ω q l
          * ((1 : Matrix (Fin (d N)) (Fin (d N)) ℝ) - m.V N * (m.V N)ᵀ) l j := rfl
  rw [h]
  exact Finset.measurable_sum _ fun l _ =>
    ((measurable_pi_apply l).comp
      ((measurable_pi_apply q).comp (m.measurable_stackEG N))).mul_const _

/-- `Y = Σ^{1/2} E⊥` is measurable. -/
theorem measurable_YhetR (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ) :
    Measurable fun ω => m.SigmaHalfR w N * m.EperpHetR N ω := by
  refine measurable_pi_lambda _ fun q => measurable_pi_lambda _ fun j => ?_
  have h : (fun ω => (m.SigmaHalfR w N * m.EperpHetR N ω) q j)
      = fun ω => ∑ l, m.SigmaHalfR w N q l * m.EperpHetR N ω l j := rfl
  rw [h]
  exact Finset.measurable_sum _ fun l _ =>
    ((measurable_pi_apply j).comp
      ((measurable_pi_apply l).comp (m.measurable_EperpHetR N))).const_mul _

/-- `Yᵀ` is measurable, in the shape `measurable_gramLamMax` composes with. The rank-one
mirror is `measurable_Yhet_transpose` (`RMT/Het/R5het.lean:404`). -/
theorem measurable_YhetR_transpose (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ)
    (N : ℕ) : Measurable fun ω => (m.SigmaHalfR w N * m.EperpHetR N ω)ᵀ := by
  refine measurable_pi_lambda _ fun j => measurable_pi_lambda _ fun q => ?_
  have h : (fun ω => (m.SigmaHalfR w N * m.EperpHetR N ω)ᵀ j q)
      = fun ω => ∑ l, m.SigmaHalfR w N q l * m.EperpHetR N ω l j := rfl
  rw [h]
  exact Finset.measurable_sum _ fun l _ =>
    ((measurable_pi_apply j).comp
      ((measurable_pi_apply l).comp (m.measurable_EperpHetR N))).const_mul _

/-- `W₀'` is measurable, through `W₀' = Y Yᵀ`. -/
theorem measurable_W0hetR (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ) :
    Measurable fun ω => m.W0hetR w N ω := by
  have h : (fun ω => m.W0hetR w N ω)
      = fun ω => (m.SigmaHalfR w N * m.EperpHetR N ω)
        * (m.SigmaHalfR w N * m.EperpHetR N ω)ᵀ := by
    funext ω
    exact m.W0hetR_eq_mul_transpose w N ω
  rw [h]
  refine measurable_pi_lambda _ fun q => measurable_pi_lambda _ fun q' => ?_
  have he : (fun ω => ((m.SigmaHalfR w N * m.EperpHetR N ω)
        * (m.SigmaHalfR w N * m.EperpHetR N ω)ᵀ) q q')
      = fun ω => ∑ l, (m.SigmaHalfR w N * m.EperpHetR N ω) q l
          * (m.SigmaHalfR w N * m.EperpHetR N ω) q' l := rfl
  rw [he]
  exact Finset.measurable_sum _ fun l _ =>
    (((measurable_pi_apply l).comp
        ((measurable_pi_apply q).comp (m.measurable_YhetR w N))).mul
      ((measurable_pi_apply l).comp
        ((measurable_pi_apply q').comp (m.measurable_YhetR w N))))

/-- `λ_max(W₀') = gramLamMax Yᵀ`, because `W₀' = Y Yᵀ = (Yᵀ)ᵀ Yᵀ`. The rank-one mirror is
`lamMax_W0het_eq_gramLamMax` (`RMT/Het/R5het.lean:415`). -/
theorem lamMax_W0hetR_eq_gramLamMax (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ)
    (N : ℕ) (ω : Ω N) :
    lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω)
      = gramLamMax (m.SigmaHalfR w N * m.EperpHetR N ω)ᵀ := by
  have h : m.W0hetR w N ω = ((m.SigmaHalfR w N * m.EperpHetR N ω)ᵀ)ᵀ
      * (m.SigmaHalfR w N * m.EperpHetR N ω)ᵀ := by
    rw [m.W0hetR_eq_mul_transpose w N ω, Matrix.transpose_transpose]
  exact lamMax_congr h _ _

/-- **The edge event is measurable.** The input of stage E3. The rank-one mirror is
`measurableSet_lamMax_W0het_le` (`RMT/Het/R5het.lean:425`). -/
theorem measurableSet_lamMax_W0hetR_le (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ)
    (N : ℕ) (t : ℝ) :
    MeasurableSet {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) ≤ t} := by
  have h : {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) ≤ t}
      = (fun ω => gramLamMax (m.SigmaHalfR w N * m.EperpHetR N ω)ᵀ) ⁻¹' Set.Iic t := by
    ext ω
    simp only [Set.mem_ofPred_eq, Set.mem_preimage, Set.mem_Iic,
      lamMax_W0hetR_eq_gramLamMax]
  rw [h]
  exact (measurable_gramLamMax.comp (m.measurable_YhetR_transpose w N)) measurableSet_Iic

end UnalignedModelR

end StackedSVD
