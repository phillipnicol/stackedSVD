/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.Prob.NoiseLaw

/-! # Stage 3, the sharp upper edge at a general law: the shared definitions

The definitions every unit of `notes/stage3_edge.md` (route C, the moment method with a
truncation) shares. The plan and the unit table are in `notes/STAGE3_CAMPAIGN.md`.

The noise matrix `Z` (entries i.i.d. of law `ν`, `NoiseLaw ν`: mean 0, variance 1, finite
fourth moment) splits as `Z = Ẑ + R + m_T J` at the level `T = d^(1/2 - a)`:
`Ẑ = truncMat ν T Z` (the entries at most `T` in absolute value, recentered), `R = discardMat T Z`
(the entries above `T`), and `m_T J = meanMat ν T n d` (the constant matrix the recentering
adds back). The law of one entry of `Ẑ` is `truncLaw ν T`, a `TruncNoiseLaw` (centered,
variance at most 1, support in `[-2T, 2T]`). The trace moment `∫ trace ((ẐᵀẐ)^k)` at
`k = momOrder C d = ⌈C log d⌉` is expanded over the closed walks of the complete bipartite
graph; `walkMult` is the multiplicity of an entry in one walk, `IsPaired` and `IsExcess` are
the two walk classes with no multiplicity 1, and `excessSum` is the weighted sum over the
excess class. `fkFactor` is the Furedi-Komlos price of one unit of vertex excess. -/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix ENNReal NNReal

namespace StackedSVD
namespace Edge

/-- The cyclic successor on `Fin k`. `Fin k` has no `OfNat 1` instance without `NeZero k`,
so the walk index is written with this map rather than with `t + 1`. -/
def cycSucc {k : ℕ} (t : Fin k) : Fin k :=
  ⟨(t.1 + 1) % k, Nat.mod_lt _ (Nat.lt_of_le_of_lt (Nat.zero_le t.1) t.2)⟩

/-- `m_T = ∫ x 1{|x| ≤ T} ∂ν`, the mean that the truncation removes. -/
noncomputable def truncMean (ν : Measure ℝ) (T : ℝ) : ℝ :=
  ∫ x, (if |x| ≤ T then x else 0) ∂ν

/-- One truncated and recentered entry. -/
noncomputable def truncMap (ν : Measure ℝ) (T : ℝ) (x : ℝ) : ℝ :=
  (if |x| ≤ T then x else 0) - truncMean ν T

/-- One discarded entry. -/
noncomputable def discardMap (T : ℝ) (x : ℝ) : ℝ := if |x| ≤ T then 0 else x

/-- The law of one truncated and recentered entry. -/
noncomputable def truncLaw (ν : Measure ℝ) (T : ℝ) : Measure ℝ := Measure.map (truncMap ν T) ν

/-- `Ẑ`, the entrywise truncation of a matrix. -/
noncomputable def truncMat (ν : Measure ℝ) (T : ℝ) {n d : ℕ}
    (Z : Matrix (Fin n) (Fin d) ℝ) : Matrix (Fin n) (Fin d) ℝ :=
  Matrix.of fun i j => truncMap ν T (Z i j)

/-- `R`, the discarded part of a matrix. -/
noncomputable def discardMat (T : ℝ) {n d : ℕ}
    (Z : Matrix (Fin n) (Fin d) ℝ) : Matrix (Fin n) (Fin d) ℝ :=
  Matrix.of fun i j => discardMap T (Z i j)

/-- `m_T J`, the constant matrix the recentering adds back. -/
noncomputable def meanMat (ν : Measure ℝ) (T : ℝ) (n d : ℕ) : Matrix (Fin n) (Fin d) ℝ :=
  Matrix.of fun _ _ => truncMean ν T

/-- The truncation level `T_N = d^(1/2 - a)`, `0 < a < 1/8` (choice 2 of the note). -/
noncomputable def truncLevel (a : ℝ) (d : ℕ) : ℝ := (d : ℝ) ^ ((1 : ℝ) / 2 - a)

/-- The truncation level is positive. -/
theorem truncLevel_pos {a : ℝ} {D : ℕ} (hD : 0 < D) : 0 < truncLevel a D :=
  Real.rpow_pos_of_pos (by exact_mod_cast hD) _

/-- The moment order `k_N = ⌈C log d⌉`. -/
noncomputable def momOrder (C : ℝ) (d : ℕ) : ℕ := ⌈C * Real.log d⌉₊

/-- The Markov constant: at `C = markovConst c ε` and `k = momOrder C d` the Markov bound
`d ((bulkEdge c + ε)/(bulkEdge c + 2ε))^k` is at most `1/d`. -/
noncomputable def markovConst (c ε : ℝ) : ℝ :=
  2 / Real.log ((bulkEdge c + 2 * ε) / (bulkEdge c + ε))

/-- The Furedi-Komlos excess factor: the price of one unit of vertex excess,
`K² (2k)¹² / min(n, d)` (Furedi and Komlos 1981, page 237, at walk length `2k`). The
exponent 12 is the one the count of unit X delivers; the assembly reads only that the factor
is at most `1/2` for `N` large (`eventually_fkFactor_le`), which holds for every polynomial
in `k`. -/
noncomputable def fkFactor (K : ℝ) (n d k : ℕ) : ℝ :=
  K ^ 2 * (2 * k : ℝ) ^ 12 / (min n d : ℝ)

/-- The noise class of the truncated law: centered, variance at most 1, bounded support.
`NoiseLaw` asks for variance exactly 1, which the truncation lowers, so units C, E and X
consume this class rather than `NoiseLaw`. -/
structure TruncNoiseLaw (ρ : Measure ℝ) (K : ℝ) : Prop where
  prob : IsProbabilityMeasure ρ
  mean : ∫ x, x ∂ρ = 0
  var_le : ∫ x, x ^ 2 ∂ρ ≤ 1
  bdd : ∀ᵐ x ∂ρ, |x| ≤ K

/-- The multiplicity of the entry `e` among the `2k` factors of the closed walk
`j 0 → i 0 → j 1 → i 1 → ... → j (k-1) → i (k-1) → j 0`. The factors are
`Y (i t) (j t)` and `Y (i t) (j (cycSucc t))`, `t : Fin k`. -/
def walkMult {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) (e : Fin n × Fin d) : ℕ :=
  ({t : Fin k | i t = e.1 ∧ j t = e.2} : Finset (Fin k)).card +
    ({t : Fin k | i t = e.1 ∧ j (cycSucc t) = e.2} : Finset (Fin k)).card

/-- Class `A`: every entry multiplicity is `0` or `2`. Weight `(∫ x² ∂ρ)^k` at any law. -/
def IsPaired {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) : Prop :=
  ∀ e : Fin n × Fin d, walkMult i j e = 0 ∨ walkMult i j e = 2

/-- Class `B`: no entry multiplicity `1`, and some entry multiplicity at least `3`. -/
def IsExcess {k n d : ℕ} (i : Fin k → Fin n) (j : Fin k → Fin d) : Prop :=
  (∀ e : Fin n × Fin d, walkMult i j e ≠ 1) ∧ ∃ e : Fin n × Fin d, 3 ≤ walkMult i j e

open scoped Classical in
/-- The weighted sum over class `B`, with the absolute moments of `ρ`. Unit X bounds it. -/
noncomputable def excessSum (ρ : Measure ℝ) (n d k : ℕ) : ℝ :=
  ∑ i : Fin k → Fin n, ∑ j : Fin k → Fin d,
    if IsExcess i j then ∏ e : Fin n × Fin d, ∫ x, |x| ^ (walkMult i j e) ∂ρ else 0

end Edge
end StackedSVD
