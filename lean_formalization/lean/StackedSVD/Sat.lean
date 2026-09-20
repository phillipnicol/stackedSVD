/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.Het.Sup
import StackedSVD.SVDStack.Weighted
import StackedSVD.RankR.Het.Sup

/-!
# Concrete Gaussian models: the hypothesis sets are satisfiable

This file answers one question an auditor must ask of any formal statement: is the hypothesis
set satisfiable? A theorem whose hypotheses no model meets is true and empty. Each theorem
below builds a model in this file and applies one of the three Gaussian endpoints to it, with
no hypothesis left over.

| Witness | Endpoint | Paper |
|---|---|---|
| `Sat.Rank1.sat_stacksvd_weighted` | `thm_stacksvd_weighted_gaussian` | `thm:stacksvd_weighted` |
| `Sat.Rank1.sat_svdstack_weighted` | `thm_svdstack_weighted_gaussian` | `thm:svdstack_weighted` |
| `Sat.RankR.sat_rank_r_stacksvd` | `thm_rank_r_stacksvd_gaussian` | `thm:rank_r_stacksvd` |

The module is in the import closure of `StackedSVD.lean`, so `scripts/check_axioms.sh` audits
these three declarations on every run and the witnesses cannot rot.

Scope. A witness settles satisfiability of one hypothesis set. It says nothing about other
models, and it changes no statement. The conclusions are not trivial at these instances: every
model here is above the detection threshold, so every limit is strictly positive.

Sources. The rank-1 construction is the vacuity probe of the D33 audit
(`notes/archive/audit_independent_D33_2026-09-02.md`); the rank-`r` one is section 5.3 of
`notes/archive/audit_independent_TrackE_wave2b_2026-09-02.md`. Both ran outside the tree; this file
brings them in.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

namespace StackedSVD

namespace Sat

/-! ## Rank one: `M` tables, `n_i N = k_i (N+1)`, `d N = l (N+1)`, so `c_i = k_i / l` -/

namespace Rank1

/-- The sample space at index `N`: one noise matrix per table. -/
abbrev Om (M : ℕ) (n : Fin M → ℕ → ℕ) (d : ℕ → ℕ) (N : ℕ) : Type :=
  (i : Fin M) → Matrix (Fin (n i N)) (Fin (d N)) ℝ

/-- Independent Gaussian noise for the `M` tables. -/
noncomputable def mu (M : ℕ) (n : Fin M → ℕ → ℕ) (d : ℕ → ℕ) (N : ℕ) : Measure (Om M n d N) :=
  Measure.pi fun i => gaussianMatrix (n i N) (d N)

instance instProbMu (M : ℕ) (n : Fin M → ℕ → ℕ) (d : ℕ → ℕ) (N : ℕ) :
    IsProbabilityMeasure (mu M n d N) := by
  unfold mu; infer_instance

/-- `n_i N = k_i (N + 1)`. -/
abbrev nk {M : ℕ} (k : Fin M → ℕ) (i : Fin M) (N : ℕ) : ℕ := k i * (N + 1)

/-- `d N = l (N + 1)`, so the aspect ratio is exactly `k_i / l` at every `N`. -/
abbrev dl (l : ℕ) (N : ℕ) : ℕ := l * (N + 1)

/-- One table: spike `θ i` on the first coordinate, Gaussian noise `ω i`. -/
noncomputable def tbl {M : ℕ} (k : Fin M → ℕ) (l : ℕ) (hk : ∀ i, 0 < k i) (hl : 0 < l)
    (θ : Fin M → ℝ) (hθ : ∀ i, 0 ≤ θ i) (i : Fin M) :
    SpikedModel (mu M (nk k) (dl l)) (nk k i) (dl l) where
  θ := θ i
  u := fun N => EuclideanSpace.single ⟨0, Nat.mul_pos (hk i) N.succ_pos⟩ 1
  v := fun N => EuclideanSpace.single ⟨0, Nat.mul_pos hl N.succ_pos⟩ 1
  Z := fun N ω => ω i
  hθ := hθ i
  hn := fun N => Nat.mul_pos (hk i) N.succ_pos
  hd := fun N => Nat.mul_pos hl N.succ_pos
  hu := fun N => by simp
  hv := fun N => by simp
  hZ := fun N => measurable_pi_apply i

/-- The `M` tables as one model with a shared `v`. -/
noncomputable def model {M : ℕ} (k : Fin M → ℕ) (l : ℕ) (hk : ∀ i, 0 < k i) (hl : 0 < l)
    (θ : Fin M → ℝ) (hθ : ∀ i, 0 ≤ θ i) :
    MultiTableModel (mu M (nk k) (dl l)) M (nk k) (dl l) where
  tbl := tbl k l hk hl θ hθ
  hv := fun _ _ _ => rfl

theorem model_joint {M : ℕ} (k : Fin M → ℕ) (l : ℕ) (hk : ∀ i, 0 < k i) (hl : 0 < l)
    (θ : Fin M → ℝ) (hθ : ∀ i, 0 ≤ θ i) :
    (model k l hk hl θ hθ).JointGaussianNoise := fun N => by
  unfold mu
  exact HasLaw.id

theorem model_θ {M : ℕ} (k : Fin M → ℕ) (l : ℕ) (hk : ∀ i, 0 < k i) (hl : 0 < l)
    (θ : Fin M → ℝ) (hθ : ∀ i, 0 ≤ θ i) (i : Fin M) :
    ((model k l hk hl θ hθ).tbl i).θ = θ i := rfl

theorem model_regime {M : ℕ} (k : Fin M → ℕ) (l : ℕ) (hk : ∀ i, 0 < k i) (hl : 0 < l)
    (θ : Fin M → ℝ) (hθ : ∀ i, 0 ≤ θ i) (c : Fin M → ℝ) (hc : ∀ i, (k i : ℝ) / l = c i)
    (i : Fin M) : ((model k l hk hl θ hθ).tbl i).Regime (c i) := by
  unfold SpikedModel.Regime
  refine ⟨?_, ?_, ?_⟩
  · refine tendsto_atTop_mono (fun N => ?_) tendsto_id
    calc N ≤ N + 1 := Nat.le_succ N
      _ = 1 * (N + 1) := (one_mul _).symm
      _ ≤ k i * (N + 1) := Nat.mul_le_mul_right _ (hk i)
  · refine tendsto_atTop_mono (fun N => ?_) tendsto_id
    calc N ≤ N + 1 := Nat.le_succ N
      _ = 1 * (N + 1) := (one_mul _).symm
      _ ≤ l * (N + 1) := Nat.mul_le_mul_right _ hl
  · have h : ∀ N : ℕ, ((nk k i N : ℕ) : ℝ) / ((dl l N : ℕ) : ℝ) = c i := by
      intro N
      rw [← hc i]
      have hN : ((N : ℝ) + 1) ≠ 0 := by positivity
      simp only [nk, dl]
      push_cast
      exact mul_div_mul_right _ _ hN
    exact tendsto_const_nhds.congr (fun N => (h N).symm)

/-! ### The instance: `M = 2`, `c = (1, 1)`, `θ = (2, 2)`, both tables above threshold -/

/-- Two tables, square, spike 2 in each: `θ_i⁴ / c_i = 16 > 1`. -/
noncomputable def mA : MultiTableModel (mu 2 (nk fun _ => 1) (dl 1)) 2 (nk fun _ => 1) (dl 1) :=
  model (fun _ => 1) 1 (fun _ => one_pos) one_pos ![2, 2] (by intro i; fin_cases i <;> norm_num)

theorem mA_θ (i : Fin 2) : (mA.tbl i).θ = ![2, 2] i := rfl

theorem regA (i : Fin 2) : (mA.tbl i).Regime ((fun _ => (1 : ℝ)) i) :=
  model_regime _ _ _ _ _ _ (fun _ => 1) (by intro i; simp) i

theorem jointA : mA.JointGaussianNoise := model_joint _ _ _ _ _ _

theorem hθA : ∃ i, (mA.tbl i).θ ≠ 0 := ⟨0, by rw [mA_θ]; norm_num⟩

/-- `thm:stacksvd_weighted` applied to `mA`: the hypothesis set of
`thm_stacksvd_weighted_gaussian` is satisfiable. -/
theorem sat_stacksvd_weighted :
    TendstoInProb (mu 2 (nk fun _ => 1) (dl 1))
      (fun N ω => mA.stackPerfW (Scalars.optWstack (fun i => (mA.tbl i).θ) (fun _ => 1)) N ω)
      (Scalars.stackSVDLimitW (fun i => (mA.tbl i).θ) (fun _ => 1)) :=
  mA.thm_stacksvd_weighted_gaussian (fun _ => 1) (fun _ => one_pos) hθA regA jointA

/-- `thm:svdstack_weighted` applied to `mA`: the hypothesis set of
`thm_svdstack_weighted_gaussian` is satisfiable. -/
theorem sat_svdstack_weighted :
    TendstoInProb (mu 2 (nk fun _ => 1) (dl 1))
      (fun N ω => mA.svdstackPerfW (optW (fun i => beta (mA.tbl i).θ 1)) N ω)
      (svdstackLimitOpt (fun i => beta (mA.tbl i).θ 1)) :=
  mA.thm_svdstack_weighted_gaussian (fun _ => 1) (fun i => beta (mA.tbl i).θ 1)
    (fun _ => one_pos) (fun _ => rfl)
    ⟨0, by rw [mA_θ]; exact beta_pos_of_thr one_pos (by norm_num)⟩ regA jointA

end Rank1

/-! ## Rank `r`: two tables, `M = 2`, `r = 2`, `θ = [[3, 3/2], [5/2, 1]]`, `R_i = 1` -/

namespace RankR

/-- `k` orthonormal columns inside `Fin D`: the first `k` columns of the identity. -/
private noncomputable def frame (D k : ℕ) (h : k ≤ D) : Matrix (Fin D) (Fin k) ℝ :=
  (1 : Matrix (Fin D) (Fin D) ℝ).submatrix id (Fin.castLE h)

private theorem frame_orth (D k : ℕ) (h : k ≤ D) : (frame D k h)ᵀ * frame D k h = 1 := by
  ext a b
  simp only [frame, Matrix.mul_apply, Matrix.transpose_apply, Matrix.submatrix_apply, id_eq,
    Matrix.one_apply]
  rw [Finset.sum_eq_single (Fin.castLE h a)]
  · by_cases hab : a = b
    · subst hab; simp
    · simp [hab]
  · intro c _ hc
    simp [hc]
  · intro hc
    exact absurd (Finset.mem_univ _) hc

/-- `n_i N = N + 3`, the same for both tables. -/
abbrev nn : Fin 2 → ℕ → ℕ := fun _ N => N + 3

/-- `d N = N + 3`, so `c_i = 1`. -/
abbrev dd : ℕ → ℕ := fun N => N + 3

/-- The sample space at index `N`. -/
abbrev Om (N : ℕ) : Type := (i : Fin 2) → Matrix (Fin (nn i N)) (Fin (dd N)) ℝ

/-- Independent Gaussian noise for the two tables. -/
noncomputable def mu (N : ℕ) : Measure (Om N) :=
  Measure.pi fun i : Fin 2 => gaussianMatrix (nn i N) (dd N)

instance instProbMu (N : ℕ) : IsProbabilityMeasure (mu N) := by
  unfold mu; infer_instance

/-- `θ = [[3, 3/2], [5/2, 1]]`: positive and strictly decreasing in each table. -/
noncomputable def th : Fin 2 → Fin 2 → ℝ := ![![3, 3/2], ![5/2, 1]]

theorem th_pos (i k : Fin 2) : 0 < th i k := by
  fin_cases i <;> fin_cases k <;> norm_num [th]

theorem th_anti (i : Fin 2) : StrictAnti (th i) := by
  intro a b hab
  fin_cases i <;> fin_cases a <;> fin_cases b <;> simp_all [th]
  norm_num

/-- Table `i` of the rank-`r` model. -/
noncomputable def tbl (i : Fin 2) : SpikedModelR mu (nn i) dd 2 where
  θ := th i
  U := fun N => frame (nn i N) 2 (by change (2:ℕ) ≤ N + 3; omega)
  V := fun N => frame (dd N) 2 (by change (2:ℕ) ≤ N + 3; omega)
  Z := fun _ ω => ω i
  hθnn := fun k => (th_pos i k).le
  hθanti := th_anti i
  hn := fun N => by change (0:ℕ) < N + 3; omega
  hd := fun N => by change (0:ℕ) < N + 3; omega
  hU := fun N => frame_orth _ _ _
  hV := fun N => frame_orth _ _ _
  hZ := fun N => measurable_pi_apply i

/-- The two tables as one exactly aligned rank-`r` model, `R_i = 1`. -/
noncomputable def mdl : UnalignedModelR mu 2 nn dd 2 (alignedRk 2 2) where
  tbl := tbl
  V := fun N => frame (dd N) 2 (by change (2:ℕ) ≤ N + 3; omega)
  R := alignedR 2 2
  hV := fun N => frame_orth _ _ _
  hR := fun i => by simp [alignedR]
  hv := fun i N => by simp [tbl, alignedR]

theorem hc : ∀ i : Fin 2, (0:ℝ) < (fun _ : Fin 2 => (1:ℝ)) i := fun _ => one_pos

theorem hR : ∀ i, mdl.R i = 1 := fun _ => rfl

theorem hreg : ∀ i, (mdl.tbl i).Regime ((fun _ : Fin 2 => (1:ℝ)) i) := by
  intro i
  refine ⟨?_, ?_, ?_⟩
  · exact tendsto_atTop_mono (fun N => Nat.le_add_right N 3) tendsto_id
  · exact tendsto_atTop_mono (fun N => Nat.le_add_right N 3) tendsto_id
  · refine tendsto_const_nhds.congr fun N => ?_
    change (1:ℝ) = ((N + 3 : ℕ) : ℝ) / ((N + 3 : ℕ) : ℝ)
    have h : ((N + 3 : ℕ) : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
    exact (div_self h).symm

theorem hG : mdl.JointGaussianNoise := fun _ => ProbabilityTheory.HasLaw.id

/-- `thm:rank_r_stacksvd` applied to `mdl`: the hypothesis set of
`thm_rank_r_stacksvd_gaussian` is satisfiable. Both components are above the threshold
(`∑_i θ_i0⁴ / c_i = 120.06`, `∑_i θ_i1⁴ / c_i = 6.06`), so both `γ_j` are strictly positive
(`0.928594873932` and `0.614906970741`, matched to `1.0e-13` against `compute_x_star` of
`theory_pred.R`). -/
theorem sat_rank_r_stacksvd :
    (∀ j : Fin 2, TendstoInProb mu
        (fun N ω => ⟪mdl.vhatStackR (fun _ => 1) j N ω, mdl.colVecG N j⟫_ℝ ^ 2)
        (Scalars.gammaR mdl.thetaAligned (fun _ => 1) j)) ∧
      TendstoInProb mu (fun N ω => mdl.frobSqStackR (fun _ => 1) N ω)
        (∑ j : Fin 2, Scalars.gammaR mdl.thetaAligned (fun _ => 1) j) :=
  mdl.thm_rank_r_stacksvd_gaussian (fun _ => 1) hc hR hreg hG (fun j => ⟨0, th_pos 0 j⟩)

end RankR

end Sat

end StackedSVD
