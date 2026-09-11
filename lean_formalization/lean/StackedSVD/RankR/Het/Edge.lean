/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.EdgeSharp
import StackedSVD.RankR.Het.Split
import StackedSVD.RankR.RMT.TableStack
import StackedSVD.RankR.StackGamma

/-!
# Stage E3: the rank-`r` heteroscedastic edge, and the tail shift of the model

Stage E3 of `notes/archive/rankr_TrackE_plan.md` (change 5 of
`notes/archive/audit_rankr_plan_E_2026-09-02.md`). Two independent parts.

## 1. The edge

`HeteroEdgeR m w c b` says `λ_max(W₀') ≤ b + ε` with probability tending to `1`, for every
`ε > 0`. It is the rank-`r` copy of `MultiTableModel.HeteroEdge` (`RMT/Het/R4het.lean:94`)
on the split block `W0hetR` of `RankR/Het/Split.lean`.

`heteroEdgeR_of_gaussian` proves it for Gaussian noise at the exact edge `MPhet.bHet c w`.
The route is the rank-one one of `MultiTableModel.heteroEdge_of_gaussian`
(`RMT/Het/EdgeSharp.lean`), with two changes only.

1. The block of `exists_block_hasLaw_hetR` has `p N` columns with `d N = p N + r`, not
   `d N - 1` columns. `tendsto_measure_lamMax_le_of_bound` was generalized from `p + 1 = d`
   to `p + r = d` for this (same stage).
2. The limit of the normalized bound needs `p N / d N → 1`, which follows from
   `p N = d N - r` and `d N → ∞`. At `r = 1` this is the rank-one `(d - 1) / d → 1`.

The scalar side (the sharp Sudakov-Fernique bound `integral_opNorm_mul_le_sharp` at
`γ⋆ = -1/s⋆`, and the identity `bd N ² / d N = γ⋆ p N / d N + ∑_i (n_i/d) w_i² γ⋆/(γ⋆ - w_i²)`)
is unchanged: it never sees `r`.

The side condition `d N = p N + r` is not a restriction on the model.
`heteroEdgeR_of_gaussian_shift` removes it with the tail shift of part 2, at the cost of
`N ↦ N + k`.

## 2. The tail shift of `UnalignedModelR`

`shift m k` reindexes the model by `N ↦ N + k`: same `θ`, same `R`, every sequence read at
`N + k`. The mirrors are `SpikedModel.shift` (`RMT/TailShift.lean:51`), `SpikedModelR.shift`
(`RankR/RMT/TableStack.lean:272`) and `RankRStack.shift` (`RankR/RMT/ShiftR.lean:55`). Every
bridge is `rfl`, so a statement proved for the shifted model transports with no rewriting of
the objects. `exists_shift_add` is the shift that makes `d (N + k) = p N + r` hold for every
`N`, which is what the block law of stage E1 asks. `HeteroLawR.of_shift` transports the
law of `RankR/StackGamma.lean:489` back from the shifted model; its two limit fields are tail
properties and its almost sure field `simpleIdxJ` enters as an argument, exactly as
`SingleTableLaw.of_shift` (`RMT/TailShift.lean:129`) handles `topSimple`.

Paper: `main_paper.tex:2294` (`eq:rank_r_model`), `:2325` (`eq:stacksvd_appXstack`).

No `sorry`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix Matrix.Norms.L2Operator InnerProductSpace ENNReal NNReal

namespace StackedSVD

namespace UnalignedModelR

variable {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} {r : ℕ} {rk : Fin M → ℕ}

/-! ### 1. The edge as a hypothesis structure -/

/-- **The heteroscedastic edge at rank `r`**, as a hypothesis structure for a general noise
law. The field says `λ_max(W₀') ≤ b + ε` with probability tending to `1`, for every `ε > 0`.
`b` is a parameter, as `c` is; `heteroEdgeR_of_gaussian` proves the structure at
`b = MPhet.bHet c w` for Gaussian noise. The rank-one mirror is `MultiTableModel.HeteroEdge`
(`RMT/Het/R4het.lean:94`). -/
structure HeteroEdgeR (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (b : ℝ) : Prop where
  /-- `λ_max(W₀') ≤ b + ε` with probability tending to `1`, for every `ε > 0`. -/
  edge : ∀ ε > 0, Tendsto (fun N => μ N
    {ω | lamMax (m.W0hetR w N ω) (m.isHermitian_W0hetR w N ω) ≤ b + ε}) atTop (𝓝 1)

/-! ### 2. Two facts the edge proof reads -/

/-- `‖Σ^{1/2}‖ ≤ √(wSqMax w) = max_i |w_i|`. `SigmaHalfR` is the diagonal matrix of the
weights, so the rank-one proof of `MultiTableModel.opNorm_SigmaHalf_le`
(`RMT/Het/R3het.lean:667`) applies verbatim; `SigmaHalfR_eq_shell` (`RankR/Het/Split.lean:109`)
is the same statement through the shell, at the cost of `[NeZero M]`. -/
theorem opNorm_SigmaHalfR_le (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ) :
    ‖m.SigmaHalfR w N‖ ≤ Real.sqrt (Scalars.wSqMax w) := by
  rw [SigmaHalfR]
  refine R3het.opNorm_diagonal_le _ (Real.sqrt_nonneg _) fun q => ?_
  rw [← Real.sqrt_sq_eq_abs]
  exact Real.sqrt_le_sqrt (Scalars.le_wSqMax w _)

/-- The stack has at least one row. The rank-one mirror is `MultiTableModel.stack_row_pos`
(`StackSVD.lean:106`). -/
theorem stack_row_pos [NeZero M] (m : UnalignedModelR μ M n d r rk) (N : ℕ) :
    0 < ∑ i, n i N := by
  have : Nonempty (Fin M) := ⟨0⟩
  exact Finset.sum_pos (fun i _ => (m.tbl i).hn N) Finset.univ_nonempty

/-! ### 3. The edge for Gaussian noise -/

/-- **The rank-`r` heteroscedastic edge for Gaussian noise** (audit change 5). For every
`ε > 0`, `λ_max(W₀') ≤ bHet c w + ε` with probability tending to `1`.

The bound sequence is the sharp Sudakov-Fernique bound of `integral_opNorm_mul_le_sharp` at
the fixed `γ⋆ = -1/s⋆`, which is admissible because `γ⋆ > max_i w_i²`
(`MPhet.sq_lt_neg_inv_sStar`):

`bd N = √(γ⋆ p N + ∑_q a_q² γ⋆/(γ⋆ - a_q²))`, `a_q = w_{blk(q)}`.

`bd N ² / d N = γ⋆ p N / d N + ∑_i (n_i/d) w_i² γ⋆/(γ⋆ - w_i²)` is an identity, so the limit
is `MPhet.edgeObjective c w γ⋆ = MPhet.bHet c w` by `p N / d N → 1` and `n_i/d → c_i`. The
rank-one mirror is `MultiTableModel.heteroEdge_of_gaussian` (`RMT/Het/EdgeSharp.lean:786`);
only the column count of the block changes, from `d N - 1` to `p N`. Paper: the bulk edge
`b(c, w)` of `thm:stacksvd_weighted` (`main_paper.tex:463`), read at rank `r` for
`thm:rank_r_stacksvd` (`:2337`) under `assum:general_noise` (`:246`). -/
theorem heteroEdgeR_of_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    {p : ℕ → ℕ} (hpd : ∀ N, d N = p N + r) :
    m.HeteroEdgeR w c (MPhet.bHet c w) := by
  classical
  set γ : ℝ := -1 / MPhet.sStar c w with hγdef
  have hγ0 : (0 : ℝ) < γ := MPhet.neg_inv_sStar_pos hc hw
  have hwγ : ∀ i, w i ^ 2 < γ := fun i => MPhet.sq_lt_neg_inv_sStar hc hw i
  -- the plumbing of stage E1
  have hd : ∀ N, 0 < d N := (m.tbl 0).hd
  choose B hB hlaw using fun N => m.exists_block_hasLaw_hetR hG N (hpd N)
  have hW : ∀ N ω, m.W0hetR w N ω = ((d N : ℝ))⁻¹ •
      ((m.SigmaHalfR w N * B N ω) * (m.SigmaHalfR w N * B N ω)ᵀ) :=
    fun N ω => m.W0hetR_eq_of_block w N ω (B N ω) (hB N ω)
  have hlawB : ∀ N, HasLaw (B N) (gaussianMatrix (∑ i, n i N) (p N)) (μ N) :=
    fun N => measurePreserving_snd.fun_comp_hasLaw (hlaw N)
  -- the stacked diagonal of the weights
  set a : (N : ℕ) → Fin (∑ i, n i N) → ℝ := fun N q => w (finSigmaFinEquiv.symm q).1 with hadef
  have hSig : ∀ N, m.SigmaHalfR w N = Matrix.diagonal (a N) := fun N => rfl
  have haγ : ∀ (N : ℕ) (q : Fin (∑ i, n i N)), (a N q) ^ 2 < γ := fun N q => hwγ _
  -- the sharp bound sequence and its square
  set bd : ℕ → ℝ := fun N => Real.sqrt (γ * ((p N : ℕ) : ℝ)
    + ∑ q : Fin (∑ i, n i N), (a N q) ^ 2 * γ / (γ - (a N q) ^ 2)) with hbddef
  have hrad : ∀ N, (0 : ℝ) ≤ γ * ((p N : ℕ) : ℝ)
      + ∑ q : Fin (∑ i, n i N), (a N q) ^ 2 * γ / (γ - (a N q) ^ 2) := by
    intro N
    refine add_nonneg (by positivity) (Finset.sum_nonneg fun q _ => ?_)
    exact div_nonneg (mul_nonneg (sq_nonneg _) hγ0.le) (sub_pos.mpr (haγ N q)).le
  have hbd0 : ∀ N, (0 : ℝ) ≤ bd N := fun N => Real.sqrt_nonneg _
  have hbdsq : ∀ N, bd N ^ 2 = γ * ((p N : ℕ) : ℝ)
      + ∑ q : Fin (∑ i, n i N), (a N q) ^ 2 * γ / (γ - (a N q) ^ 2) := by
    intro N
    simp only [hbddef]
    exact Real.sq_sqrt (hrad N)
  -- the expectation bound, once `p N` is positive
  have hd2 : ∀ᶠ N in atTop, r + 1 ≤ d N :=
    (hreg 0).2.1.eventually (eventually_ge_atTop (r + 1))
  have hbd : ∀ᶠ N in atTop, (∫ B', ‖m.SigmaHalfR w N * B'‖
      ∂(gaussianMatrix (∑ i, n i N) (p N))) ≤ bd N := by
    filter_upwards [hd2] with N hN
    have hpp : 0 < p N := by have := hpd N; omega
    have h := integral_opNorm_mul_le_sharp (m.stack_row_pos N) hpp (a N) hγ0 (haγ N)
    rw [hSig N]
    simpa only [hbddef] using h
  -- regrouping the diagonal sum by table
  have hregroup : ∀ N, ∑ q : Fin (∑ i, n i N), (a N q) ^ 2 * γ / (γ - (a N q) ^ 2)
      = ∑ i, (n i N : ℝ) * (w i ^ 2 * γ / (γ - w i ^ 2)) := by
    intro N
    simp only [hadef]
    rw [← Equiv.sum_comp finSigmaFinEquiv]
    simp only [Equiv.symm_apply_apply]
    rw [Fintype.sum_sigma]
    simp [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  -- the limit of the normalized square: an exact identity, then `p/d → 1` and `n_i/d → c_i`
  have hval : γ + ∑ i, c i * (w i ^ 2 * γ / (γ - w i ^ 2)) = MPhet.bHet c w := by
    rw [← MPhet.edgeObjective_sStar_eq_bHet hc hw, ← hγdef]
    have hterm : ∀ i, c i * (w i ^ 2 * γ / (γ - w i ^ 2))
        = γ * (c i * w i ^ 2 / (γ - w i ^ 2)) := fun i => by ring
    simp only [hterm, MPhet.edgeObjective]
    rw [mul_add, mul_one, Finset.mul_sum]
  -- `p N / d N = 1 - r / d N → 1`; at `r = 1` this is the rank-one `(d - 1)/d → 1`
  have hpc : ∀ N, ((p N : ℕ) : ℝ) = (d N : ℝ) - (r : ℝ) := fun N => by
    rw [hpd N]; push_cast; ring
  have hone : Tendsto (fun N => (((p N : ℕ) : ℝ)) / (d N : ℝ)) atTop (𝓝 1) := by
    have hdR : Tendsto (fun N => (d N : ℝ)) atTop atTop :=
      tendsto_natCast_atTop_atTop.comp (hreg 0).2.1
    have h : Tendsto (fun N => 1 - (r : ℝ) * ((d N : ℝ))⁻¹) atTop (𝓝 1) := by
      have h0 := (tendsto_inv_atTop_zero.comp hdR).const_mul (r : ℝ)
      simpa using tendsto_const_nhds.sub h0
    refine h.congr fun N => ?_
    have hdN : ((d N : ℝ)) ≠ 0 := by
      have hpos : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast hd N
      exact hpos.ne'
    rw [hpc N]
    field_simp
  have hlim : Tendsto (fun N => bd N ^ 2 / (d N : ℝ)) atTop (𝓝 (MPhet.bHet c w)) := by
    rw [← hval]
    have h1 : Tendsto (fun N => γ * ((((p N : ℕ) : ℝ)) / (d N : ℝ))) atTop (𝓝 γ) := by
      simpa using hone.const_mul γ
    have h2 : Tendsto (fun N => ∑ i, ((n i N : ℝ) / (d N : ℝ)) * (w i ^ 2 * γ / (γ - w i ^ 2)))
        atTop (𝓝 (∑ i, c i * (w i ^ 2 * γ / (γ - w i ^ 2)))) :=
      tendsto_finsetSum _ fun i _ => (hreg i).2.2.mul_const _
    refine (h1.add h2).congr fun N => ?_
    rw [hbdsq N, hregroup N, add_div, Finset.sum_div]
    congr 1
    · ring
    · exact Finset.sum_congr rfl fun i _ => by ring
  refine ⟨fun ε hε => ?_⟩
  exact tendsto_measure_lamMax_le_of_bound μ (fun N => m.SigmaHalfR w N)
    (Real.sqrt_nonneg _) (MPhet.bHet_pos hc hw).le (fun N => m.opNorm_SigmaHalfR_le w N)
    bd hbd0 hbd hlim B (m.W0hetR w) hW (m.isHermitian_W0hetR w) hlawB m.stack_row_pos hd
    (fun N => (hpd N).symm) (hreg 0).2.1 hε

/-! ### 4. The tail shift of the model -/

/-- **The shifted model.** `UnalignedModelR` reindexed by `N ↦ N + k`: the same `θ` in every
table, the same alignment matrices `R_i`, every sequence read at `N + k`. The probability
spaces move with the index, so `Ω` becomes `fun N => Ω (N + k)` and `μ` becomes
`fun N => μ (N + k)`. The mirrors are `SpikedModel.shift` (`RMT/TailShift.lean:51`),
`SpikedModelR.shift` (`RankR/RMT/TableStack.lean:272`) and `RankRStack.shift`
(`RankR/RMT/ShiftR.lean:55`). -/
def shift (m : UnalignedModelR μ M n d r rk) (k : ℕ) :
    UnalignedModelR (fun N => μ (N + k)) M (fun i N => n i (N + k)) (fun N => d (N + k)) r
      rk where
  tbl := fun i => (m.tbl i).shift k
  V := fun N => m.V (N + k)
  R := m.R
  hV := fun N => m.hV (N + k)
  hR := m.hR
  hv := fun i N => m.hv i (N + k)

section ShiftBridges

variable (m : UnalignedModelR μ M n d r rk) (k N : ℕ)

@[simp] theorem shift_tbl (i : Fin M) : (m.shift k).tbl i = (m.tbl i).shift k := rfl

@[simp] theorem shift_V : (m.shift k).V N = m.V (N + k) := rfl

@[simp] theorem shift_R : (m.shift k).R = m.R := rfl

@[simp] theorem shift_colVecG (l : Fin r) :
    (m.shift k).colVecG N l = m.colVecG (N + k) l := rfl

@[simp] theorem shift_stackZG (ω : Ω (N + k)) :
    (m.shift k).stackZG N ω = m.stackZG (N + k) ω := rfl

@[simp] theorem shift_stackEG (ω : Ω (N + k)) :
    (m.shift k).stackEG N ω = m.stackEG (N + k) ω := rfl

@[simp] theorem shift_stackXG (ω : Ω (N + k)) :
    (m.shift k).stackXG N ω = m.stackXG (N + k) ω := rfl

@[simp] theorem shift_signalFactorG :
    (m.shift k).signalFactorG N = m.signalFactorG (N + k) := rfl

@[simp] theorem shift_stackXW (w : Fin M → ℝ) (ω : Ω (N + k)) :
    (m.shift k).stackXW w N ω = m.stackXW w (N + k) ω := rfl

@[simp] theorem shift_stackGramW (w : Fin M → ℝ) (ω : Ω (N + k)) :
    (m.shift k).stackGramW w N ω = m.stackGramW w (N + k) ω := rfl

@[simp] theorem shift_SigmaHalfR (w : Fin M → ℝ) :
    (m.shift k).SigmaHalfR w N = m.SigmaHalfR w (N + k) := rfl

@[simp] theorem shift_EperpHetR (ω : Ω (N + k)) :
    (m.shift k).EperpHetR N ω = m.EperpHetR (N + k) ω := rfl

@[simp] theorem shift_W0hetR (w : Fin M → ℝ) (ω : Ω (N + k)) :
    (m.shift k).W0hetR w N ω = m.W0hetR w (N + k) ω := rfl

end ShiftBridges

section ShiftBridgesAligned

variable (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ) (k N : ℕ)

@[simp] theorem shift_thetaAligned : (m.shift k).thetaAligned = m.thetaAligned := rfl

@[simp] theorem shift_stackXJ (j : Fin r) (ω : Ω (N + k)) :
    (m.shift k).stackXJ c j N ω = m.stackXJ c j (N + k) ω := rfl

@[simp] theorem shift_stackGramJ (j : Fin r) (ω : Ω (N + k)) :
    (m.shift k).stackGramJ c j N ω = m.stackGramJ c j (N + k) ω := rfl

@[simp] theorem shift_vhatStackR (j : Fin r) (ω : Ω (N + k)) :
    (m.shift k).vhatStackR c j N ω = m.vhatStackR c j (N + k) ω := rfl

@[simp] theorem shift_stackOverlapJ (j : Fin r) (ω : Ω (N + k)) :
    (m.shift k).stackOverlapJ c j N ω = m.stackOverlapJ c j (N + k) ω := rfl

end ShiftBridgesAligned

/-! ### 5. Transport of the hypotheses -/

/-- The joint Gaussian noise law is stated index by index, so it shifts by evaluation. The
mirror is `SpikedModelR.shift_GaussianNoise` (`RankR/RMT/TableStack.lean:317`). -/
theorem shift_JointGaussianNoise (m : UnalignedModelR μ M n d r rk) (k : ℕ)
    (hG : m.JointGaussianNoise) : (m.shift k).JointGaussianNoise := fun N => hG (N + k)

/-- The regime of one table is a tail property. -/
theorem shift_Regime (m : UnalignedModelR μ M n d r rk) (k : ℕ) (i : Fin M) {c : ℝ}
    (h : (m.tbl i).Regime c) : ((m.shift k).tbl i).Regime c :=
  (m.tbl i).shift_Regime k h

/-- Transport of the `IsProbabilityMeasure` instance to the shifted spaces. The mirror is
`SpikedModel.isProbabilityMeasure_shift` (`RMT/TailShift.lean:100`). -/
theorem isProbabilityMeasure_shift [∀ N, IsProbabilityMeasure (μ N)] (k : ℕ) :
    ∀ N, IsProbabilityMeasure (μ (N + k)) := fun _ => inferInstance

/-- `Tendsto d atTop atTop` gives a shift after which `d (N + k) = p N + r` for every `N`,
with `0 < p N`. This is the shape the block law of stage E1
(`exists_block_hasLaw_hetR`, `RankR/Het/Split.lean:263`) and `heteroEdgeR_of_gaussian` take.
The mirror is `RankRStack.exists_shift_lt` (`RankR/RMT/ShiftR.lean:98`). The model `_m` is
not read; it is a binder so that a consumer can write `m.exists_shift_add`. -/
theorem exists_shift_add (_m : UnalignedModelR μ M n d r rk) (hdtop : Tendsto d atTop atTop) :
    ∃ k, ∃ p : ℕ → ℕ, (∀ N, d (N + k) = p N + r) ∧ ∀ N, 0 < p N := by
  obtain ⟨k, hk⟩ := eventually_atTop.mp (hdtop.eventually_ge_atTop (r + 1))
  refine ⟨k, fun N => d (N + k) - r, fun N => ?_, fun N => ?_⟩ <;>
    · have h : r + 1 ≤ d (N + k) := hk (N + k) (Nat.le_add_left k N)
      dsimp only
      omega

/-! ### 6. The edge with no side condition -/

/-- The edge is a tail property: it holds for the model as soon as it holds for a shift of
it. Every object of the statement is `rfl`-equal to its shifted twin
(`shift_W0hetR`), and the two `IsHermitian` proofs match by proof irrelevance, so only
`SpikedModel.tendsto_of_shift` (`RMT/TailShift.lean:110`) is needed. -/
theorem HeteroEdgeR.of_shift {m : UnalignedModelR μ M n d r rk} {w c : Fin M → ℝ} {b : ℝ}
    (k : ℕ) (h : (m.shift k).HeteroEdgeR w c b) : m.HeteroEdgeR w c b :=
  ⟨fun ε hε => SpikedModel.tendsto_of_shift k (h.edge ε hε)⟩

/-- **The edge for Gaussian noise, through the shift** (the form stage E9 asks for). The
side condition `d N = p N + r` of `heteroEdgeR_of_gaussian` is met after a shift, because
`d N → ∞`. Paper: the bulk edge `b(c, w)` of `thm:stacksvd_weighted` (`main_paper.tex:463`)
at rank `r`. -/
theorem heteroEdgeR_of_gaussian_shift [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    ∃ k, (m.shift k).HeteroEdgeR w c (MPhet.bHet c w) := by
  obtain ⟨k, p, hp, -⟩ := m.exists_shift_add (hreg 0).2.1
  have : ∀ N, IsProbabilityMeasure (μ (N + k)) := isProbabilityMeasure_shift k
  exact ⟨k, heteroEdgeR_of_gaussian (m.shift k) w c hc hw
    (fun i => m.shift_Regime k i (hreg i)) (m.shift_JointGaussianNoise k hG) (p := p) hp⟩

/-- **The edge for Gaussian noise, with no side condition.** The `hpd` of
`heteroEdgeR_of_gaussian` and the shift of `heteroEdgeR_of_gaussian_shift` both disappear:
the shift is applied and then removed by `HeteroEdgeR.of_shift`. This is the statement a
consumer should cite. Paper: the bulk edge `b(c, w)` of `thm:stacksvd_weighted`
(`main_paper.tex:463`) at rank `r`. -/
theorem heteroEdgeR_of_gaussian_tail [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    m.HeteroEdgeR w c (MPhet.bHet c w) := by
  obtain ⟨k, hk⟩ := m.heteroEdgeR_of_gaussian_shift w c hc hw hreg hG
  exact HeteroEdgeR.of_shift k hk

/-! ### 7. The law transports back from the shift -/

/-- **`HeteroLawR` is a tail property, except for `simpleIdxJ`.** The two limit fields
`align` and `crossProj` come from the shifted model through
`SpikedModel.tendstoInProb_of_shift` (`RMT/TailShift.lean:116`); every object in them is
`rfl`-equal to its shifted twin (section 4). `simpleIdxJ` is a statement about one `N` at a
time, so it enters as the argument `hsimple`. The mirror is `SingleTableLaw.of_shift`
(`RMT/TailShift.lean:129`). The field names are those of `RankR/StackGamma.lean:489` after
the Track E audit (`align`, `crossProj`, `simpleIdxJ`). -/
theorem HeteroLawR.of_shift {m : UnalignedModelR μ M n d r (alignedRk M r)} {c : Fin M → ℝ}
    (k : ℕ) (law : (m.shift k).HeteroLawR c)
    (hsimple : ∀ (j : Fin r) (N : ℕ), ∀ᵐ ω ∂(μ N), SimpleIdx (m.stackGramJ c j N ω)
      (m.isHermitian_stackGramJ c j N ω) (Scalars.ellR m.thetaAligned c j)) :
    m.HeteroLawR c := by
  refine ⟨fun j => ?_, fun j l hl => ?_, hsimple⟩
  · exact SpikedModel.tendstoInProb_of_shift k (law.align j)
  · exact SpikedModel.tendstoInProb_of_shift k (law.crossProj j l hl)

end UnalignedModelR

end StackedSVD
