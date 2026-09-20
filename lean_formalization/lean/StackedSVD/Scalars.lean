/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.StackSVD
import StackedSVD.SVDStack.Defs

/-!
# Class (A): the scalar results

Milestone M2 of `PLAN.md`. Every result here is an inequality or an identity between finitely
many real numbers. No matrix and no probability appears. The review note is
`notes/archive/class_A.md`.

## Content

1. `svdTerm`, the summand of `S`; the `betaSq` facts that link it to `beta`.
2. `Sval_beta` and `Sval_beta_nonneg`, in terms of `svdTerm`. `Sval` and
   `svdstackLimitOpt = S/(S+1)`, the limit of optimally weighted svdstack
   (`thm:svdstack_weighted`), are defined once in `SVDStack/Defs.lean`; this file and
   `SVDStack/Weighted.lean` both read them from there
   (`notes/archive/audit_weighted_2026-08-30.md`, section 4).
3. `binaryStackSVDLimit S θ c`, the limit of binary-weighted stacksvd on a subset (`cor.2`),
   and `binaryStackSVDLimit_eq_stackSVDLimit`, the identity that says it is the unweighted
   stacksvd limit of the sub-collection.
4. `gW` and `stackSVDLimitW`, the limit of optimally weighted stacksvd
   (`thm:stacksvd_weighted`): the unique root in `(0,1)` of `gW = 1`, and `0` below the
   threshold `∑ θ_i⁴/c_i > 1`.
5. `svdstackOpt_le_binary` (`thm:stacksvd_binary_optimal_svd_stack`).
6. `svdstackOpt_le_stackSVDLimitW` and `stackSVDLimit_le_stackSVDLimitW` (`prop:dominance`).

## Conventions

`θ` enters only through `θ ^ 2` and `θ ^ 4`, so no result needs `0 ≤ θ i`; the numeric scan
of `notes/archive/class_A.md` confirms this. Every result needs `0 < c i`, since `c i = 0` makes
`βᵢ² = 1` and `gW 0` infinite.
-/

open Finset

namespace StackedSVD
namespace Scalars

variable {M : ℕ}

/-! ### Two elementary inequalities

The first is the paper's `eq:helpfulIneq`. The second is the shape in which the two
"dominates" theorems use it: `s ↦ s/(s+1)` is increasing, and the target value is written as
`A/(A+B)`.
-/

/-- `eq:helpfulIneq`: `x/z ≤ (x+y)/(y+z)` for `0 < y`, `0 < z` and `x ≤ z`. -/
theorem div_le_add_div_add {x y z : ℝ} (hy : 0 < y) (hz : 0 < z) (hxz : x ≤ z) :
    x / z ≤ (x + y) / (y + z) := by
  rw [div_le_div_iff₀ hz (by linarith)]
  nlinarith

/-- `s/(s+1) ≤ A/(A+B)` whenever `s * B ≤ A`. Both limits of the paper have this shape:
`s` is `S` of `thm:svdstack_weighted` and `A/(A+B)` is the stacksvd limit. -/
theorem div_succ_le_div_add (s A B : ℝ) (hs : 0 ≤ s) (hA : 0 < A) (hB : 0 < B)
    (h : s * B ≤ A) : s / (s + 1) ≤ A / (A + B) := by
  rw [div_le_div_iff₀ (by linarith) (by linarith)]
  nlinarith

/-- The strict form of `div_succ_le_div_add`. -/
theorem div_succ_lt_div_add (s A B : ℝ) (hs : 0 ≤ s) (hA : 0 < A) (hB : 0 < B)
    (h : s * B < A) : s / (s + 1) < A / (A + B) := by
  rw [div_lt_div_iff₀ (by linarith) (by linarith)]
  nlinarith

/-! ### `betaSq` and the summand of `S` -/

/-- The summand of `S` of `thm:svdstack_weighted` written in `(θ, c)`:
`β²/(1-β²) = (θ⁴-c)/(θ²+c)` above the threshold, `0` below it. -/
noncomputable def svdTerm (θ c : ℝ) : ℝ := if c < θ ^ 4 then (θ ^ 4 - c) / (θ ^ 2 + c) else 0

theorem svdTerm_nonneg {θ c : ℝ} (hc : 0 < c) : 0 ≤ svdTerm θ c := by
  unfold svdTerm
  split_ifs with h
  · exact div_nonneg (by linarith) (by nlinarith [sq_nonneg θ, hc])
  · exact le_rfl

/-- `betaSq` is nonnegative for every `θ` and `c`: above the threshold the numerator and
the denominator are both nonnegative, and below it the value is `0`. -/
theorem betaSq_nonneg (θ c : ℝ) : 0 ≤ betaSq θ c := by
  unfold betaSq
  split_ifs with h
  · exact div_nonneg (by linarith) (by positivity)
  · exact le_rfl

theorem one_sub_betaSq_pos {θ c : ℝ} (hc : 0 < c) : 0 < 1 - betaSq θ c := by
  unfold betaSq
  split_ifs with h
  · have hd : (0:ℝ) < θ ^ 4 + θ ^ 2 := by nlinarith [sq_nonneg θ]
    have : (θ ^ 4 - c) / (θ ^ 4 + θ ^ 2) < 1 := by
      rw [div_lt_one hd]
      nlinarith [sq_nonneg θ]
    linarith
  · norm_num

/-- `β²/(1-β²) = svdTerm θ c`. This is the identity the paper uses to write `S` in `(θ, c)`. -/
theorem betaSq_div_one_sub {θ c : ℝ} (hc : 0 < c) :
    betaSq θ c / (1 - betaSq θ c) = svdTerm θ c := by
  unfold betaSq svdTerm
  split_ifs with h
  · have hd : (0:ℝ) < θ ^ 4 + θ ^ 2 := by nlinarith [sq_nonneg θ]
    have h2 : (0:ℝ) < θ ^ 2 + c := by nlinarith [sq_nonneg θ]
    have key : (1:ℝ) - (θ ^ 4 - c) / (θ ^ 4 + θ ^ 2) = (θ ^ 2 + c) / (θ ^ 4 + θ ^ 2) := by
      rw [eq_div_iff hd.ne', sub_mul, div_mul_cancel₀ _ hd.ne']
      ring
    rw [key, div_eq_div_iff (div_ne_zero h2.ne' hd.ne') h2.ne']
    ring
  · norm_num

/-- `β = √(β²)` squares back, so `beta` and `betaSq` agree under a square. -/
theorem beta_sq (θ c : ℝ) : beta θ c ^ 2 = betaSq θ c :=
  Real.sq_sqrt (betaSq_nonneg θ c)

/-! ### `S` and the optimally weighted svdstack limit (`thm:svdstack_weighted`) -/

/-- `S` in the paper's data `(θ, c)`. -/
theorem Sval_beta (θ c : Fin M → ℝ) (hc : ∀ i, 0 < c i) :
    Sval (fun i => beta (θ i) (c i)) = ∑ i, svdTerm (θ i) (c i) := by
  unfold Sval
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [beta_sq, betaSq_div_one_sub (hc i)]

theorem Sval_beta_nonneg (θ c : Fin M → ℝ) (hc : ∀ i, 0 < c i) :
    0 ≤ Sval (fun i => beta (θ i) (c i)) := by
  rw [Sval_beta θ c hc]
  exact Finset.sum_nonneg fun i _ => svdTerm_nonneg (hc i)

/-! ### Binary-weighted stacksvd (`cor.2`) -/

/-- The limit of stacksvd run on the sub-collection `S ⊆ [M]` (`cor.2`). The guard is the
detectability condition of `prop:stacksvd_general` for the sub-collection. -/
noncomputable def binaryStackSVDLimit (S : Finset (Fin M)) (θ c : Fin M → ℝ) : ℝ :=
  if (∑ i ∈ S, θ i ^ 2) ^ 2 > ∑ i ∈ S, c i then
    ((∑ i ∈ S, θ i ^ 2) ^ 2 - ∑ i ∈ S, c i) / ((∑ i ∈ S, θ i ^ 2) ^ 2 + ∑ i ∈ S, θ i ^ 2)
  else 0

/-- A sum over `S` read through any enumeration `Fin k ≃ S`. -/
theorem sum_equiv_subtype {k : ℕ} {S : Finset (Fin M)} (e : Fin k ≃ {x // x ∈ S})
    (f : Fin M → ℝ) : ∑ j : Fin k, f (e j) = ∑ i ∈ S, f i := by
  rw [← Finset.sum_coe_sort S f]
  exact Fintype.sum_equiv e _ _ fun _ => rfl

/-- `cor.2` is `prop:stacksvd_general` for the sub-collection: the binary-weights value on `S`
is `stackSVDLimit` of `θ` and `c` restricted to `S`, read through any enumeration of `S`. -/
theorem binaryStackSVDLimit_eq_stackSVDLimit {k : ℕ} (S : Finset (Fin M)) (θ c : Fin M → ℝ)
    (e : Fin k ≃ {x // x ∈ S}) :
    binaryStackSVDLimit S θ c = stackSVDLimit (fun j => θ (e j)) (fun j => c (e j)) := by
  have hθ : ∑ j : Fin k, θ (e j) ^ 2 = ∑ i ∈ S, θ i ^ 2 :=
    sum_equiv_subtype e fun i => θ i ^ 2
  have hc : ∑ j : Fin k, c (e j) = ∑ i ∈ S, c i := sum_equiv_subtype e c
  simp only [binaryStackSVDLimit, stackSVDLimit, hθ, hc]
  split_ifs with h
  · rw [show (∑ i ∈ S, θ i ^ 2) * ((∑ i ∈ S, θ i ^ 2) + 1)
        = (∑ i ∈ S, θ i ^ 2) ^ 2 + ∑ i ∈ S, θ i ^ 2 by ring]
  · rfl

/-- Above the threshold the sum of the `c_i` is below the square of the sum of the `θ_i²`:
`∑_S c_i < ∑_S θ_i⁴ = ∑_S (θ_i²)² ≤ (∑_S θ_i²)²`. So the guard of `binaryStackSVDLimit`
holds as soon as every kept table is above its own threshold. -/
theorem sum_lt_sq_sum_of_above {S : Finset (Fin M)} {θ c : Fin M → ℝ}
    (hS : ∀ i ∈ S, c i < θ i ^ 4) (hne : S.Nonempty) :
    ∑ i ∈ S, c i < (∑ i ∈ S, θ i ^ 2) ^ 2 := by
  have h1 : ∑ i ∈ S, c i < ∑ i ∈ S, θ i ^ 4 := Finset.sum_lt_sum_of_nonempty hne hS
  have h2 : ∑ i ∈ S, (θ i ^ 2) ^ 2 ≤ (∑ i ∈ S, θ i ^ 2) ^ 2 :=
    Finset.sum_sq_le_sq_sum_of_nonneg fun i _ => sq_nonneg _
  have h3 : ∑ i ∈ S, θ i ^ 4 = ∑ i ∈ S, (θ i ^ 2) ^ 2 :=
    Finset.sum_congr rfl fun i _ => by ring
  linarith [h3 ▸ h1]

/-- `cor.2`, the closed form of binary-weighted stacksvd. The `if` of
`binaryStackSVDLimit` resolves because every kept table is above its own threshold. -/
theorem binaryStackSVDLimit_eq {S : Finset (Fin M)} {θ c : Fin M → ℝ}
    (hS : ∀ i ∈ S, c i < θ i ^ 4) (hne : S.Nonempty) :
    binaryStackSVDLimit S θ c =
      ((∑ i ∈ S, θ i ^ 2) ^ 2 - ∑ i ∈ S, c i) / ((∑ i ∈ S, θ i ^ 2) ^ 2 + ∑ i ∈ S, θ i ^ 2) := by
  unfold binaryStackSVDLimit
  exact if_pos (sum_lt_sq_sum_of_above hS hne)

/-! ### `thm:stacksvd_binary_optimal_svd_stack`

The paper splits on `∑ c_i` versus `1` in three cases. One argument covers all three:
`(∑_j y_j)(∑_i c_i(1-c_i)/y_i) ≥ ∑_i c_i(1-c_i) = ∑ c_i - ∑ c_i²` and `∑ c_i² ≤ (∑ c_i)²`,
with the second strict as soon as two tables are kept.
-/

/-- `∑ a_i ≤ (∑ y_j)(∑ a_i/y_i)` for nonnegative `a` and positive `y`. -/
theorem sum_le_sum_mul_sum_div {S : Finset (Fin M)} {a y : Fin M → ℝ}
    (ha : ∀ i ∈ S, 0 ≤ a i) (hy : ∀ i ∈ S, 0 < y i) :
    ∑ i ∈ S, a i ≤ (∑ j ∈ S, y j) * ∑ i ∈ S, a i / y i := by
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum fun i hi => ?_
  have hyi : y i ≠ 0 := (hy i hi).ne'
  have h1 : y i ≤ ∑ j ∈ S, y j := Finset.single_le_sum (fun j hj => (hy j hj).le) hi
  have h2 : 0 ≤ a i / y i := div_nonneg (ha i hi) (hy i hi).le
  calc a i = y i * (a i / y i) := by field_simp
    _ ≤ (∑ j ∈ S, y j) * (a i / y i) := mul_le_mul_of_nonneg_right h1 h2

/-- Step 1 of the paper's proof: `(θ⁴-c)/(θ²+c) = (θ²-c) - c(1-c)/(θ²+c)`. -/
theorem svdTerm_decomp {θ c : ℝ} (hc : 0 < c) (h : c < θ ^ 4) :
    svdTerm θ c = (θ ^ 2 - c) - c * (1 - c) / (θ ^ 2 + c) := by
  have hd : (0:ℝ) < θ ^ 2 + c := by nlinarith [sq_nonneg θ]
  rw [svdTerm, if_pos h]
  field_simp
  ring

/-- `∑ c_i² < (∑ c_i)²` for two or more positive `c_i`. -/
theorem sum_sq_lt_sq_sum {S : Finset (Fin M)} {c : Fin M → ℝ} (hc : ∀ i ∈ S, 0 < c i)
    (hcard : 2 ≤ S.card) : ∑ i ∈ S, c i ^ 2 < (∑ i ∈ S, c i) ^ 2 := by
  obtain ⟨i0, hi0⟩ := Finset.card_pos.mp (by omega : 0 < S.card)
  have hS'ne : (S.erase i0).Nonempty := by
    rw [← Finset.card_pos, Finset.card_erase_of_mem hi0]
    omega
  have hsplit : c i0 + ∑ i ∈ S.erase i0, c i = ∑ i ∈ S, c i := Finset.add_sum_erase S c hi0
  have hsplit2 : c i0 ^ 2 + ∑ i ∈ S.erase i0, c i ^ 2 = ∑ i ∈ S, c i ^ 2 :=
    Finset.add_sum_erase S (fun i => c i ^ 2) hi0
  have hb : 0 < ∑ i ∈ S.erase i0, c i :=
    Finset.sum_pos (fun i hi => hc i (Finset.mem_of_mem_erase hi)) hS'ne
  have hle : ∑ i ∈ S.erase i0, c i ^ 2 ≤ (∑ i ∈ S.erase i0, c i) ^ 2 :=
    Finset.sum_sq_le_sq_sum_of_nonneg fun i hi => (hc i (Finset.mem_of_mem_erase hi)).le
  have hc0 : 0 < c i0 := hc i0 hi0
  rw [← hsplit, ← hsplit2]
  nlinarith

/-- The core inequality of `thm:stacksvd_binary_optimal_svd_stack`, in the form the
comparison `div_succ_le_div_add` consumes: `S · (T + C) ≤ (T² - C) - ((∑ c)² - ∑ c²)`
with `T = ∑_S θ_i²` and `C = ∑_S c_i`. The last bracket is nonnegative, and positive as soon
as two tables are kept. -/
theorem binary_core {S : Finset (Fin M)} {θ c : Fin M → ℝ} (hc : ∀ i ∈ S, 0 < c i)
    (hc1 : ∀ i ∈ S, c i ≤ 1) (hS : ∀ i ∈ S, c i < θ i ^ 4) :
    (∑ i ∈ S, svdTerm (θ i) (c i)) * ((∑ i ∈ S, θ i ^ 2) + ∑ i ∈ S, c i)
      ≤ ((∑ i ∈ S, θ i ^ 2) ^ 2 - ∑ i ∈ S, c i)
        - ((∑ i ∈ S, c i) ^ 2 - ∑ i ∈ S, c i ^ 2) := by
  have hy : ∀ i ∈ S, 0 < θ i ^ 2 + c i := fun i hi => by nlinarith [sq_nonneg (θ i), hc i hi]
  have hdec : ∑ i ∈ S, svdTerm (θ i) (c i)
      = ((∑ i ∈ S, θ i ^ 2) - ∑ i ∈ S, c i)
        - ∑ i ∈ S, c i * (1 - c i) / (θ i ^ 2 + c i) := by
    rw [show ((∑ i ∈ S, θ i ^ 2) - ∑ i ∈ S, c i) = ∑ i ∈ S, (θ i ^ 2 - c i) from
        (Finset.sum_sub_distrib _ _).symm, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun i hi => svdTerm_decomp (hc i hi) (hS i hi)
  have hsum_a : ∑ i ∈ S, c i * (1 - c i) = (∑ i ∈ S, c i) - ∑ i ∈ S, c i ^ 2 := by
    rw [show ((∑ i ∈ S, c i) - ∑ i ∈ S, c i ^ 2) = ∑ i ∈ S, (c i - c i ^ 2) from
        (Finset.sum_sub_distrib _ _).symm]
    exact Finset.sum_congr rfl fun i _ => by ring
  have hsum_y : ∑ i ∈ S, (θ i ^ 2 + c i) = (∑ i ∈ S, θ i ^ 2) + ∑ i ∈ S, c i :=
    Finset.sum_add_distrib
  have hkey := sum_le_sum_mul_sum_div (S := S) (a := fun i => c i * (1 - c i))
    (y := fun i => θ i ^ 2 + c i)
    (fun i hi => mul_nonneg (hc i hi).le (by linarith [hc1 i hi])) hy
  rw [hsum_a, hsum_y] at hkey
  rw [hdec]
  nlinarith [hkey]

/-- `thm:stacksvd_binary_optimal_svd_stack`: with every kept table at `c_i ≤ 1`,
binary-weighted stacksvd is at least optimally weighted svdstack. `S` is the set of tables
that the binary weighting keeps, `w_i = 1{β_i > 0}`. -/
theorem svdstackOpt_le_binary {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hc1 : ∀ i, c i < θ i ^ 4 → c i ≤ 1) :
    svdstackLimitOpt (fun i => beta (θ i) (c i))
      ≤ binaryStackSVDLimit (Finset.univ.filter fun i => c i < θ i ^ 4) θ c := by
  classical
  set S : Finset (Fin M) := Finset.univ.filter fun i => c i < θ i ^ 4 with hSdef
  have hmem : ∀ i, i ∈ S ↔ c i < θ i ^ 4 := by
    intro i; rw [hSdef, Finset.mem_filter]; simp
  have hSval : Sval (fun i => beta (θ i) (c i)) = ∑ i ∈ S, svdTerm (θ i) (c i) := by
    rw [Sval_beta θ c hc]
    refine (Finset.sum_subset (Finset.subset_univ S) ?_).symm
    intro i _ hi
    rw [svdTerm, if_neg (fun h => hi ((hmem i).mpr h))]
  have hs0 : 0 ≤ Sval (fun i => beta (θ i) (c i)) := Sval_beta_nonneg θ c hc
  rcases S.eq_empty_or_nonempty with hemp | hne
  · rw [svdstackLimitOpt, hSval, hemp]
    simp [binaryStackSVDLimit]
  · have hSpos : ∀ i ∈ S, 0 < c i := fun i _ => hc i
    have hSabove : ∀ i ∈ S, c i < θ i ^ 4 := fun i hi => (hmem i).mp hi
    have hSc1 : ∀ i ∈ S, c i ≤ 1 := fun i hi => hc1 i (hSabove i hi)
    have hguard : ∑ i ∈ S, c i < (∑ i ∈ S, θ i ^ 2) ^ 2 := sum_lt_sq_sum_of_above hSabove hne
    have hC : 0 < ∑ i ∈ S, c i := Finset.sum_pos hSpos hne
    have hT : 0 ≤ ∑ i ∈ S, θ i ^ 2 := Finset.sum_nonneg fun i _ => sq_nonneg _
    have hsq : ∑ i ∈ S, c i ^ 2 ≤ (∑ i ∈ S, c i) ^ 2 :=
      Finset.sum_sq_le_sq_sum_of_nonneg fun i hi => (hSpos i hi).le
    have hcore := binary_core hSpos hSc1 hSabove
    rw [binaryStackSVDLimit_eq hSabove hne, svdstackLimitOpt, hSval]
    rw [show (∑ i ∈ S, θ i ^ 2) ^ 2 + ∑ i ∈ S, θ i ^ 2
        = ((∑ i ∈ S, θ i ^ 2) ^ 2 - ∑ i ∈ S, c i)
          + ((∑ i ∈ S, θ i ^ 2) + ∑ i ∈ S, c i) by ring]
    refine div_succ_le_div_add _ _ _ (by rw [← hSval]; exact hs0) (by linarith) (by linarith) ?_
    linarith

/-- The strict half of `thm:stacksvd_binary_optimal_svd_stack`: two kept tables give a strict
improvement. In the paper this is the hypothesis `β₂ > 0`. -/
theorem svdstackOpt_lt_binary {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hc1 : ∀ i, c i < θ i ^ 4 → c i ≤ 1)
    (hcard : 2 ≤ (Finset.univ.filter fun i => c i < θ i ^ 4).card) :
    svdstackLimitOpt (fun i => beta (θ i) (c i))
      < binaryStackSVDLimit (Finset.univ.filter fun i => c i < θ i ^ 4) θ c := by
  classical
  set S : Finset (Fin M) := Finset.univ.filter fun i => c i < θ i ^ 4 with hSdef
  have hmem : ∀ i, i ∈ S ↔ c i < θ i ^ 4 := by
    intro i; rw [hSdef, Finset.mem_filter]; simp
  have hne : S.Nonempty := Finset.card_pos.mp (by omega)
  have hSval : Sval (fun i => beta (θ i) (c i)) = ∑ i ∈ S, svdTerm (θ i) (c i) := by
    rw [Sval_beta θ c hc]
    refine (Finset.sum_subset (Finset.subset_univ S) ?_).symm
    intro i _ hi
    rw [svdTerm, if_neg (fun h => hi ((hmem i).mpr h))]
  have hs0 : 0 ≤ Sval (fun i => beta (θ i) (c i)) := Sval_beta_nonneg θ c hc
  have hSpos : ∀ i ∈ S, 0 < c i := fun i _ => hc i
  have hSabove : ∀ i ∈ S, c i < θ i ^ 4 := fun i hi => (hmem i).mp hi
  have hSc1 : ∀ i ∈ S, c i ≤ 1 := fun i hi => hc1 i (hSabove i hi)
  have hguard : ∑ i ∈ S, c i < (∑ i ∈ S, θ i ^ 2) ^ 2 := sum_lt_sq_sum_of_above hSabove hne
  have hC : 0 < ∑ i ∈ S, c i := Finset.sum_pos hSpos hne
  have hT : 0 ≤ ∑ i ∈ S, θ i ^ 2 := Finset.sum_nonneg fun i _ => sq_nonneg _
  have hsq : ∑ i ∈ S, c i ^ 2 < (∑ i ∈ S, c i) ^ 2 := sum_sq_lt_sq_sum hSpos hcard
  have hcore := binary_core hSpos hSc1 hSabove
  rw [binaryStackSVDLimit_eq hSabove hne, svdstackLimitOpt, hSval]
  rw [show (∑ i ∈ S, θ i ^ 2) ^ 2 + ∑ i ∈ S, θ i ^ 2
      = ((∑ i ∈ S, θ i ^ 2) ^ 2 - ∑ i ∈ S, c i)
        + ((∑ i ∈ S, θ i ^ 2) + ∑ i ∈ S, c i) by ring]
  refine div_succ_lt_div_add _ _ _ (by rw [← hSval]; exact hs0) (by linarith) (by linarith) ?_
  linarith

/-! ### Optimally weighted stacksvd (`thm:stacksvd_weighted`)

`gW θ c x = ∑ θ_i⁴(1-x)/(c_i + xθ_i²)` is the paper's `f(x) + 1`. It is continuous and
strictly decreasing on `[0,1]`, equals `∑ θ_i⁴/c_i` at `0` and `0` at `1`, so it meets `1`
exactly once in `(0,1)` when `∑ θ_i⁴/c_i > 1` (`eq:assumption4`), and never when
`∑ θ_i⁴/c_i ≤ 1`.
-/

/-- `f(x) + 1` of `prop:dominance`; the limit of optimally weighted stacksvd is the root of
`gW θ c x = 1`. -/
noncomputable def gW (θ c : Fin M → ℝ) (x : ℝ) : ℝ :=
  ∑ i, θ i ^ 4 * (1 - x) / (c i + x * θ i ^ 2)

theorem gW_denom_pos {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i) {x : ℝ} (hx : 0 ≤ x) (i : Fin M) :
    0 < c i + x * θ i ^ 2 :=
  add_pos_of_pos_of_nonneg (hc i) (mul_nonneg hx (sq_nonneg _))

theorem gW_one (θ c : Fin M → ℝ) : gW θ c 1 = 0 := by simp [gW]

theorem gW_zero (θ c : Fin M → ℝ) : gW θ c 0 = ∑ i, θ i ^ 4 / c i := by simp [gW]

theorem gW_nonneg {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i) {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x ≤ 1) :
    0 ≤ gW θ c x :=
  Finset.sum_nonneg fun i _ =>
    div_nonneg (mul_nonneg (by positivity) (by linarith)) (gW_denom_pos hc hx0 i).le

theorem gW_continuousOn {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i) :
    ContinuousOn (gW θ c) (Set.Icc 0 1) := by
  unfold gW
  refine continuousOn_finsetSum _ fun i _ => ContinuousOn.div ?_ ?_ ?_
  · exact (continuous_const.mul (continuous_const.sub continuous_id)).continuousOn
  · exact (continuous_const.add (continuous_id.mul continuous_const)).continuousOn
  · exact fun x hx => (gW_denom_pos hc hx.1 i).ne'

theorem gW_strictAntiOn {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i) (hθ : ∃ i, θ i ≠ 0) :
    StrictAntiOn (gW θ c) (Set.Icc 0 1) := by
  rintro x hx y hy hxy
  obtain ⟨i0, hi0⟩ := hθ
  refine Finset.sum_lt_sum (fun i _ => ?_) ⟨i0, Finset.mem_univ _, ?_⟩
  · rw [div_le_div_iff₀ (gW_denom_pos hc hy.1 i) (gW_denom_pos hc hx.1 i)]
    nlinarith [mul_nonneg (by positivity : (0:ℝ) ≤ θ i ^ 4)
      (mul_nonneg (by linarith : (0:ℝ) ≤ y - x)
        (add_pos_of_pos_of_nonneg (hc i) (sq_nonneg (θ i))).le)]
  · rw [div_lt_div_iff₀ (gW_denom_pos hc hy.1 i0) (gW_denom_pos hc hx.1 i0)]
    have h4 : (0:ℝ) < θ i0 ^ 4 := by positivity
    nlinarith [mul_pos h4 (mul_pos (by linarith : (0:ℝ) < y - x)
      (add_pos_of_pos_of_nonneg (hc i0) (sq_nonneg (θ i0))))]

theorem exists_ne_zero_of_thr {θ c : Fin M → ℝ} (hthr : 1 < ∑ i, θ i ^ 4 / c i) :
    ∃ i, θ i ≠ 0 := by
  by_contra h
  push Not at h
  have : ∑ i, θ i ^ 4 / c i = 0 :=
    Finset.sum_eq_zero fun i _ => by rw [h i]; norm_num
  linarith

/-- Existence and uniqueness of the root of `gW θ c x = 1` in `(0,1)`, under the paper's
detectability condition `∑ θ_i⁴/c_i > 1` (`eq:assumption4`). -/
theorem existsUnique_root {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hthr : 1 < ∑ i, θ i ^ 4 / c i) :
    ∃! x : ℝ, x ∈ Set.Ioo (0:ℝ) 1 ∧ gW θ c x = 1 := by
  have hanti := gW_strictAntiOn hc (exists_ne_zero_of_thr hthr)
  have hmem : (1:ℝ) ∈ Set.Icc (gW θ c 1) (gW θ c 0) := by
    rw [gW_one, gW_zero]
    exact ⟨by norm_num, hthr.le⟩
  obtain ⟨x, hx, hxv⟩ :=
    intermediate_value_Icc' (by norm_num : (0:ℝ) ≤ 1) (gW_continuousOn hc) hmem
  have hx0 : x ≠ 0 := by
    intro h
    rw [h, gW_zero] at hxv
    linarith
  have hx1 : x ≠ 1 := by
    intro h
    rw [h, gW_one] at hxv
    norm_num at hxv
  refine ⟨x, ⟨⟨lt_of_le_of_ne hx.1 (Ne.symm hx0), lt_of_le_of_ne hx.2 hx1⟩, hxv⟩, ?_⟩
  rintro y ⟨hy, hyv⟩
  exact hanti.injOn (Set.Ioo_subset_Icc_self hy) hx (by rw [hyv, hxv])

open Classical in
/-- The limit of optimally weighted stacksvd (`thm:stacksvd_weighted`): the unique root of
`gW θ c x = 1` in `(0,1)` when there is one, and `0` otherwise. Below the threshold
`∑ θ_i⁴/c_i > 1` there is no root, so the definition returns the paper's `0`
(`stackSVDLimitW_eq_zero`). -/
noncomputable def stackSVDLimitW (θ c : Fin M → ℝ) : ℝ :=
  if h : ∃! x : ℝ, x ∈ Set.Ioo (0:ℝ) 1 ∧ gW θ c x = 1 then h.choose else 0

theorem stackSVDLimitW_nonneg (θ c : Fin M → ℝ) : 0 ≤ stackSVDLimitW θ c := by
  unfold stackSVDLimitW
  split_ifs with h
  · exact h.choose_spec.1.1.1.le
  · exact le_rfl

theorem stackSVDLimitW_spec {θ c : Fin M → ℝ}
    (h : ∃! x : ℝ, x ∈ Set.Ioo (0 : ℝ) 1 ∧ gW θ c x = 1) :
    stackSVDLimitW θ c ∈ Set.Ioo (0:ℝ) 1 ∧ gW θ c (stackSVDLimitW θ c) = 1 := by
  unfold stackSVDLimitW
  split_ifs
  exact h.choose_spec.1

/-- Below the detectability threshold the limit is `0`. -/
theorem stackSVDLimitW_eq_zero {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (h : ∑ i, θ i ^ 4 / c i ≤ 1) : stackSVDLimitW θ c = 0 := by
  unfold stackSVDLimitW
  rw [dif_neg]
  rintro ⟨x, ⟨hx, hxv⟩, -⟩
  by_cases hθ : ∃ i, θ i ≠ 0
  · have hanti := gW_strictAntiOn hc hθ
    have := hanti (Set.left_mem_Icc.mpr (by norm_num)) (Set.Ioo_subset_Icc_self hx) hx.1
    rw [gW_zero] at this
    linarith
  · push Not at hθ
    have : gW θ c x = 0 := Finset.sum_eq_zero fun i _ => by rw [hθ i]; norm_num
    rw [this] at hxv
    norm_num at hxv

/-- The comparison lemma that both halves of `prop:dominance` use: `gW` is strictly
decreasing and equals `1` at the limit, so any `y` with `gW θ c y ≥ 1` is below the limit. -/
theorem le_stackSVDLimitW {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hthr : 1 < ∑ i, θ i ^ 4 / c i) {y : ℝ} (hy : y ∈ Set.Icc (0 : ℝ) 1)
    (hgy : 1 ≤ gW θ c y) : y ≤ stackSVDLimitW θ c := by
  obtain ⟨hx, hgx⟩ := stackSVDLimitW_spec (existsUnique_root hc hthr)
  by_contra hlt
  push Not at hlt
  have := gW_strictAntiOn hc (exists_ne_zero_of_thr hthr)
    (Set.Ioo_subset_Icc_self hx) hy hlt
  linarith

/-! ### `prop:dominance` -/

/-- The value of one summand of `gW` at `x = s/(s+1)`. -/
theorem gW_term_at_div {θ c s : ℝ} (hc : 0 < c) (hs : 0 ≤ s) :
    θ ^ 4 * (1 - s / (s + 1)) / (c + s / (s + 1) * θ ^ 2)
      = θ ^ 4 / ((s + 1) * c + s * θ ^ 2) := by
  have h1 : (0:ℝ) < s + 1 := by linarith
  have hd : (0:ℝ) < (s + 1) * c + s * θ ^ 2 :=
    add_pos_of_pos_of_nonneg (mul_pos h1 hc) (mul_nonneg hs (sq_nonneg _))
  have hden : (0:ℝ) < c + s / (s + 1) * θ ^ 2 :=
    add_pos_of_pos_of_nonneg hc (mul_nonneg (div_nonneg hs h1.le) (sq_nonneg _))
  rw [div_eq_div_iff hden.ne' hd.ne']
  field_simp
  ring

/-- The termwise inequality `eq:dominance_termwise_ineq` of the paper. -/
theorem dominance_term {θ c s : ℝ} (hc : 0 < c) (hs : 0 < s) (hle : svdTerm θ c ≤ s) :
    svdTerm θ c / s ≤ θ ^ 4 * (1 - s / (s + 1)) / (c + s / (s + 1) * θ ^ 2) := by
  have h1 : (0:ℝ) < s + 1 := by linarith
  have hd : (0:ℝ) < (s + 1) * c + s * θ ^ 2 :=
    add_pos_of_pos_of_nonneg (mul_pos h1 hc) (mul_nonneg hs.le (sq_nonneg _))
  rw [gW_term_at_div hc hs.le]
  unfold svdTerm at hle ⊢
  split_ifs with h
  · have hy : (0:ℝ) < θ ^ 2 + c := by nlinarith [sq_nonneg θ]
    rw [if_pos h] at hle
    have hle' : θ ^ 4 - c ≤ s * (θ ^ 2 + c) := by
      rw [div_le_iff₀ hy] at hle
      linarith
    rw [div_div, div_le_div_iff₀ (by positivity) hd]
    nlinarith [mul_nonneg hc.le (by linarith : (0:ℝ) ≤ s * (θ ^ 2 + c) - (θ ^ 4 - c))]
  · rw [zero_div]
    positivity

/-- `prop:dominance`, first half: optimally weighted stacksvd is at least optimally weighted
svdstack. No threshold hypothesis is needed: below the svdstack threshold `S = 0` and the
left side is `0`. -/
theorem svdstackOpt_le_stackSVDLimitW {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i) :
    svdstackLimitOpt (fun i => beta (θ i) (c i)) ≤ stackSVDLimitW θ c := by
  set s := Sval (fun i => beta (θ i) (c i)) with hsdef
  have hSval : s = ∑ i, svdTerm (θ i) (c i) := by rw [hsdef, Sval_beta θ c hc]
  have hs0 : 0 ≤ s := Sval_beta_nonneg θ c hc
  rcases eq_or_lt_of_le hs0 with hz | hpos
  · rw [svdstackLimitOpt, ← hsdef, ← hz]
    simpa using stackSVDLimitW_nonneg θ c
  · -- some table is above its own threshold, so the weighted threshold holds
    have hex : ∃ i ∈ Finset.univ, (0:ℝ) < svdTerm (θ i) (c i) := by
      by_contra hcon
      push Not at hcon
      have : ∑ i, svdTerm (θ i) (c i) ≤ 0 :=
        Finset.sum_nonpos fun i hi => hcon i hi
      rw [← hSval] at this
      linarith
    obtain ⟨i0, -, hi0⟩ := hex
    have habove : c i0 < θ i0 ^ 4 := by
      by_contra h
      rw [svdTerm, if_neg h] at hi0
      exact lt_irrefl 0 hi0
    have hthr : 1 < ∑ i, θ i ^ 4 / c i := by
      have h1 : 1 < θ i0 ^ 4 / c i0 := (one_lt_div (hc i0)).mpr habove
      have h2 : θ i0 ^ 4 / c i0 ≤ ∑ i, θ i ^ 4 / c i :=
        Finset.single_le_sum (f := fun i => θ i ^ 4 / c i)
          (fun i _ => div_nonneg (by positivity) (hc i).le) (Finset.mem_univ i0)
      linarith
    have hmem : s / (s + 1) ∈ Set.Icc (0:ℝ) 1 := by
      constructor
      · positivity
      · rw [div_le_one (by linarith)]
        linarith
    have hgoal : 1 ≤ gW θ c (s / (s + 1)) := by
      have hterm : ∀ i, svdTerm (θ i) (c i) / s
          ≤ θ i ^ 4 * (1 - s / (s + 1)) / (c i + s / (s + 1) * θ i ^ 2) := by
        intro i
        refine dominance_term (hc i) hpos ?_
        rw [hSval]
        exact Finset.single_le_sum (fun j _ => svdTerm_nonneg (hc j)) (Finset.mem_univ i)
      calc (1:ℝ) = (∑ i, svdTerm (θ i) (c i)) / s := by rw [← hSval]; field_simp
        _ = ∑ i, svdTerm (θ i) (c i) / s := Finset.sum_div _ _ _
        _ ≤ gW θ c (s / (s + 1)) := Finset.sum_le_sum fun i _ => hterm i
    exact le_trans (le_of_eq rfl) (le_stackSVDLimitW hc hthr hmem hgoal)

/-- `prop:dominance`, second half: optimally weighted stacksvd is at least unweighted
stacksvd. The step is Cauchy-Schwarz in Engel form applied to
`∑ (θ_i²)²/(c_i(T²+T) + θ_i²(T²-C))`. -/
theorem stackSVDLimit_le_stackSVDLimitW {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i) :
    stackSVDLimit θ c ≤ stackSVDLimitW θ c := by
  unfold stackSVDLimit
  split_ifs with h
  · set T := ∑ i, θ i ^ 2 with hTdef
    set C := ∑ i, c i with hCdef
    have hM : 0 < M := by
      rcases Nat.eq_zero_or_pos M with rfl | h'
      · rw [hTdef, hCdef] at h; simp at h
      · exact h'
    have hne : (Finset.univ : Finset (Fin M)).Nonempty :=
      Finset.univ_nonempty_iff.mpr ⟨⟨0, hM⟩⟩
    have hC : 0 < C := Finset.sum_pos (fun i _ => hc i) hne
    have hT0 : 0 ≤ T := Finset.sum_nonneg fun i _ => sq_nonneg _
    have hT : 0 < T := by nlinarith
    have hTT : 0 < T ^ 2 + T := by positivity
    have hTTne : (T ^ 2 + T) ≠ 0 := hTT.ne'
    have hpow : ∀ x : ℝ, (x ^ 2) ^ 2 = x ^ 4 := fun x => by ring
    -- the weighted threshold, by Cauchy-Schwarz
    have hcs0 : T ^ 2 / C ≤ ∑ i, θ i ^ 4 / c i := by
      have h0 := Finset.sq_sum_div_le_sum_sq_div (Finset.univ : Finset (Fin M))
        (fun i => θ i ^ 2) (g := c) (fun i _ => hc i)
      simp only [hpow] at h0
      rw [← hTdef, ← hCdef] at h0
      exact h0
    have hthr : 1 < ∑ i, θ i ^ 4 / c i := by
      have h1 : 1 < T ^ 2 / C := (one_lt_div hC).mpr h
      linarith
    -- the unweighted limit as a point of `[0,1]`
    set u := (T ^ 2 - C) / (T ^ 2 + T) with hudef
    have hu : u ∈ Set.Icc (0:ℝ) 1 := by
      constructor
      · rw [hudef]
        apply div_nonneg (by linarith) hTT.le
      · rw [hudef, div_le_one hTT]
        linarith
    -- `gW θ c u ≥ 1`, by Cauchy-Schwarz in Engel form
    have hbpos : ∀ i ∈ (Finset.univ : Finset (Fin M)),
        0 < c i * (T ^ 2 + T) + θ i ^ 2 * (T ^ 2 - C) := fun i _ =>
      add_pos_of_pos_of_nonneg (mul_pos (hc i) hTT)
        (mul_nonneg (sq_nonneg _) (by linarith))
    have hsumb : ∑ i, (c i * (T ^ 2 + T) + θ i ^ 2 * (T ^ 2 - C)) = T ^ 2 * (T + C) := by
      rw [Finset.sum_add_distrib, ← Finset.sum_mul, ← Finset.sum_mul, ← hTdef, ← hCdef]
      ring
    have hterm : ∀ i : Fin M, θ i ^ 4 * (1 - u) / (c i + u * θ i ^ 2)
        = (T + C) * (θ i ^ 2) ^ 2 / (c i * (T ^ 2 + T) + θ i ^ 2 * (T ^ 2 - C)) := by
      intro i
      have hden : 0 < c i + u * θ i ^ 2 :=
        add_pos_of_pos_of_nonneg (hc i) (mul_nonneg hu.1 (sq_nonneg _))
      have hbi : 0 < c i * (T ^ 2 + T) + θ i ^ 2 * (T ^ 2 - C) := hbpos i (Finset.mem_univ i)
      rw [div_eq_div_iff hden.ne' hbi.ne', hudef]
      field_simp
      ring
    have hgu : 1 ≤ gW θ c u := by
      have hcs := Finset.sq_sum_div_le_sum_sq_div (Finset.univ : Finset (Fin M))
        (fun i => θ i ^ 2) (g := fun i => c i * (T ^ 2 + T) + θ i ^ 2 * (T ^ 2 - C)) hbpos
      rw [hsumb, ← hTdef] at hcs
      have hTC : 0 < T + C := by linarith
      have hfrac : T ^ 2 / (T ^ 2 * (T + C)) = 1 / (T + C) := by
        rw [mul_comm]
        field_simp
      rw [hfrac] at hcs
      have hgeq : gW θ c u
          = (T + C) * ∑ i, (θ i ^ 2) ^ 2 / (c i * (T ^ 2 + T) + θ i ^ 2 * (T ^ 2 - C)) := by
        unfold gW
        rw [Finset.mul_sum]
        refine Finset.sum_congr rfl fun i _ => ?_
        rw [hterm i, mul_div_assoc]
      rw [hgeq]
      calc (1:ℝ) = (T + C) * (1 / (T + C)) := by field_simp
        _ ≤ (T + C) * ∑ i, (θ i ^ 2) ^ 2 / (c i * (T ^ 2 + T) + θ i ^ 2 * (T ^ 2 - C)) :=
            mul_le_mul_of_nonneg_left hcs hTC.le
    have := le_stackSVDLimitW hc hthr hu hgu
    rw [hudef] at this
    rw [show T * (T + 1) = T ^ 2 + T by ring]
    exact this
  · exact stackSVDLimitW_nonneg θ c

/-! ### `thm:simple_thm1` (`cor.1`): equal `θ`, equal `c`

`main_paper.tex:321`. Both halves are stated in the paper's own closed form. The stacksvd
half specializes `stackSVDLimit`. The svdstack half needs the spectrum of the equicorrelated
`A_β = β₀² 11ᵀ + (1-β₀²) I`, which this section computes: `λ_max = 1 + (M-1)β₀²` and
`(βᵀ v_max)² = M β₀²`, so `svdstackLimit = M β₀²/(1 + (M-1)β₀²)`.
-/

section EqualTables

open scoped InnerProductSpace Matrix

variable {M : ℕ}

/-! The first four helpers of this section are public since the 2026-09-02 dedupe: `Remarks.lean`
kept primed copies of them. -/

theorem euclid_norm_sq (x : EuclideanSpace ℝ (Fin M)) : ‖x‖ ^ 2 = ∑ i, x i ^ 2 := by
  rw [EuclideanSpace.norm_sq_eq]
  exact Finset.sum_congr rfl fun i _ => by rw [Real.norm_eq_abs, sq_abs]

theorem toOp_apply (A : Matrix (Fin M) (Fin M) ℝ) (x : EuclideanSpace ℝ (Fin M))
    (j : Fin M) : toOp A x j = ∑ k, A j k * x k := rfl

/-- `lamMax` is the Rayleigh quotient at `vMax`, so a uniform bound on the quadratic form at
unit vectors bounds it. -/
theorem lamMax_le_of_quadForm {d : ℕ} (hd : 0 < d) (A : Matrix (Fin d) (Fin d) ℝ)
    (hA : A.IsHermitian) {L : ℝ}
    (h : ∀ x : EuclideanSpace ℝ (Fin d), ‖x‖ = 1 → ⟪toOp A x, x⟫_ℝ ≤ L) :
    lamMax A hA ≤ L := by
  have hmem := mem_topSpace_vMax hd A hA
  have hnv : ‖vMax A hA‖ = 1 := norm_vMax hd A hA
  have hray : ⟪toOp A (vMax A hA), vMax A hA⟫_ℝ = lamMax A hA := by
    rw [toOp_of_mem_topSpace hmem, real_inner_smul_left, real_inner_self_eq_norm_sq, hnv]
    ring
  rw [← hray]
  exact h _ hnv

/-- The quadratic form at a unit vector is a lower bound for `lamMax`. -/
theorem le_lamMax_of_unit {d : ℕ} (A : Matrix (Fin d) (Fin d) ℝ) (hA : A.IsHermitian)
    (x : EuclideanSpace ℝ (Fin d)) (hx : ‖x‖ = 1) : ⟪toOp A x, x⟫_ℝ ≤ lamMax A hA := by
  have h := inner_toOp_self_le A hA x
  rwa [hx, one_pow, mul_one] at h

/-- One row of `A_β x` for the equicorrelated `A_β`. -/
private theorem Abeta_const_row (β₀ : ℝ) (x : EuclideanSpace ℝ (Fin M)) (j : Fin M) :
    toOp (Abeta (fun _ : Fin M => β₀)) x j = β₀ ^ 2 * (∑ k, x k) + (1 - β₀ ^ 2) * x j := by
  rw [toOp_apply]
  have hterm : ∀ k : Fin M, Abeta (fun _ : Fin M => β₀) j k * x k
      = β₀ ^ 2 * x k + (if j = k then (1 - β₀ ^ 2) * x k else 0) := by
    intro k
    have hA : Abeta (fun _ : Fin M => β₀) j k
        = β₀ * β₀ + (if j = k then 1 - β₀ ^ 2 else 0) := by
      simp [Abeta, Matrix.vecMulVec_apply, Matrix.diagonal_apply]
    rw [hA]
    split_ifs with h <;> ring
  rw [Finset.sum_congr rfl fun k _ => hterm k, Finset.sum_add_distrib, ← Finset.mul_sum,
    Finset.sum_ite_eq _ j (fun k => (1 - β₀ ^ 2) * x k), if_pos (Finset.mem_univ j)]

/-- The quadratic form of the equicorrelated `A_β`. -/
private theorem quadForm_Abeta_const (β₀ : ℝ) (x : EuclideanSpace ℝ (Fin M)) :
    ⟪toOp (Abeta (fun _ : Fin M => β₀)) x, x⟫_ℝ
      = β₀ ^ 2 * (∑ i, x i) ^ 2 + (1 - β₀ ^ 2) * ∑ i, x i ^ 2 := by
  rw [real_inner_eq_dotProduct]
  have hd : (WithLp.ofLp (toOp (Abeta (fun _ : Fin M => β₀)) x)) ⬝ᵥ WithLp.ofLp x
      = ∑ i, (β₀ ^ 2 * (∑ k, x k) + (1 - β₀ ^ 2) * x i) * x i :=
    Finset.sum_congr rfl fun i _ => by rw [← Abeta_const_row β₀ x i]
  rw [hd]
  have hsplit : ∀ i : Fin M, (β₀ ^ 2 * (∑ k, x k) + (1 - β₀ ^ 2) * x i) * x i
      = (β₀ ^ 2 * (∑ k, x k)) * x i + (1 - β₀ ^ 2) * x i ^ 2 := fun i => by ring
  rw [Finset.sum_congr rfl fun i _ => hsplit i, Finset.sum_add_distrib, ← Finset.mul_sum,
    ← Finset.mul_sum]
  ring

/-- `λ_max(A_β) = 1 + (M-1)β₀²` for the equicorrelated `A_β`. The quadratic form at a unit `x`
is `β₀²(∑ x_i)² + (1-β₀²)`, at most `1 + (M-1)β₀²` by `(∑ x_i)² ≤ M ∑ x_i²`, and equal to it
at the uniform vector `x = 1/√M`. -/
theorem lamMax_Abeta_const (hM : 0 < M) (β₀ : ℝ) :
    lamMax (Abeta (fun _ : Fin M => β₀)) (isHermitian_Abeta _)
      = 1 + ((M : ℝ) - 1) * β₀ ^ 2 := by
  have hMpos : (0 : ℝ) < M := by exact_mod_cast hM
  have hspos : 0 < Real.sqrt M := Real.sqrt_pos.mpr hMpos
  have hsq : Real.sqrt M * Real.sqrt M = M := Real.mul_self_sqrt hMpos.le
  have hexp : 1 + ((M : ℝ) - 1) * β₀ ^ 2 = β₀ ^ 2 * (M : ℝ) + (1 - β₀ ^ 2) * 1 := by ring
  refine le_antisymm (lamMax_le_of_quadForm hM _ _ fun x hx => ?_) ?_
  · rw [quadForm_Abeta_const]
    have hxsq : ∑ i, x i ^ 2 = 1 := by
      have h := euclid_norm_sq x
      rw [hx, one_pow] at h
      exact h.symm
    have hcs : (∑ i, x i) ^ 2 ≤ (M : ℝ) := by
      have h := sq_sum_le_card_mul_sum_sq (s := (Finset.univ : Finset (Fin M)))
        (f := fun i => x i)
      rw [Finset.card_univ, Fintype.card_fin, hxsq, mul_one] at h
      exact h
    have hkey : β₀ ^ 2 * (∑ i, x i) ^ 2 ≤ β₀ ^ 2 * (M : ℝ) :=
      mul_le_mul_of_nonneg_left hcs (sq_nonneg β₀)
    rw [hxsq, hexp]
    linarith
  · have hMne : (M : ℝ) ≠ 0 := hMpos.ne'
    have hr2 : ((Real.sqrt M)⁻¹) ^ 2 = ((M : ℝ))⁻¹ := by
      rw [inv_pow, Real.sq_sqrt hMpos.le]
    set u : EuclideanSpace ℝ (Fin M) :=
      WithLp.toLp 2 (fun _ : Fin M => (Real.sqrt M)⁻¹) with hudef
    have huapp : ∀ i : Fin M, u i = (Real.sqrt M)⁻¹ := fun i => by simp [hudef]
    have husum2 : (∑ i, u i) ^ 2 = (M : ℝ) := by
      simp only [huapp]
      rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_pow, hr2]
      field_simp
    have husq : ∑ i, u i ^ 2 = 1 := by
      simp only [huapp]
      rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, hr2]
      field_simp
    have hunorm : ‖u‖ = 1 := by
      have h := euclid_norm_sq u
      rw [husq] at h
      rw [← Real.sqrt_sq (norm_nonneg u), h, Real.sqrt_one]
    have hle := le_lamMax_of_unit _ (isHermitian_Abeta (fun _ : Fin M => β₀)) u hunorm
    rw [quadForm_Abeta_const, husum2, husq] at hle
    rw [hexp]
    linarith

/-- `(βᵀ v_max(A_β))² = M β₀²` for the equicorrelated `A_β`. For `β₀ ≠ 0` the eigenvalue
equation forces every coordinate of the top eigenvector to be `(∑_k v_k)/M`, so `(∑_k v_k)²
= M`. For `β₀ = 0` both sides are `0`. -/
private theorem inner_beta_vMax_const (hM : 0 < M) (β₀ : ℝ) :
    ((fun _ : Fin M => β₀) ⬝ᵥ
        WithLp.ofLp (vMax (Abeta (fun _ : Fin M => β₀)) (isHermitian_Abeta _))) ^ 2
      = (M : ℝ) * β₀ ^ 2 := by
  have hMpos : (0 : ℝ) < M := by exact_mod_cast hM
  rcases eq_or_ne β₀ 0 with rfl | hβ
  · simp
  · set v := vMax (Abeta (fun _ : Fin M => β₀)) (isHermitian_Abeta _) with hvdef
    have hmem : v ∈ topSpace (Abeta (fun _ : Fin M => β₀)) (isHermitian_Abeta _) :=
      mem_topSpace_vMax hM _ _
    have hnv : ‖v‖ = 1 := norm_vMax hM _ _
    have heig := toOp_of_mem_topSpace hmem
    have hlam := lamMax_Abeta_const hM β₀
    have hrow : ∀ j : Fin M, β₀ ^ 2 * (∑ k, v k) + (1 - β₀ ^ 2) * v j
        = (1 + ((M : ℝ) - 1) * β₀ ^ 2) * v j := by
      intro j
      rw [← Abeta_const_row β₀ v j, heig, hlam]
      simp
    have hsumeq : ∀ j : Fin M, (∑ k, v k) = (M : ℝ) * v j := by
      intro j
      have h0 : β₀ ^ 2 * ((∑ k, v k) - (M : ℝ) * v j) = 0 := by linear_combination hrow j
      rcases mul_eq_zero.mp h0 with h1 | h1
      · exact absurd h1 (pow_ne_zero 2 hβ)
      · linarith
    have hvj : ∀ j : Fin M, v j = (∑ k, v k) / M := by
      intro j
      rw [hsumeq j]
      field_simp
    have hnorm1 : ∑ j, v j ^ 2 = 1 := by
      have h := euclid_norm_sq v
      rw [hnv, one_pow] at h
      exact h.symm
    have hS : (∑ k, v k) ^ 2 = (M : ℝ) := by
      have hsum2 : ∑ j : Fin M, v j ^ 2 = (∑ k, v k) ^ 2 / M := by
        have hterm : ∀ j : Fin M, v j ^ 2 = (∑ k, v k) ^ 2 / (M : ℝ) ^ 2 := by
          intro j
          rw [hvj j, div_pow]
        rw [Finset.sum_congr rfl fun j _ => hterm j, Finset.sum_const, Finset.card_univ,
          Fintype.card_fin, nsmul_eq_mul]
        field_simp
      rw [hsum2] at hnorm1
      field_simp at hnorm1
      exact hnorm1
    have hdot : ((fun _ : Fin M => β₀) ⬝ᵥ WithLp.ofLp v) = β₀ * ∑ k, v k := by
      change ∑ i, β₀ * (WithLp.ofLp v) i = β₀ * ∑ k, v k
      rw [← Finset.mul_sum]
    rw [hdot, mul_pow, hS]
    ring

/-- The svdstack limit for equal `β_i = β₀`: `M β₀²/(1 + (M-1)β₀²)`. -/
theorem svdstackLimit_const (hM : 0 < M) (β₀ : ℝ) :
    svdstackLimit (fun _ : Fin M => β₀)
      = (M : ℝ) * β₀ ^ 2 / (1 + ((M : ℝ) - 1) * β₀ ^ 2) := by
  unfold svdstackLimit
  rw [inner_beta_vMax_const hM β₀, lamMax_Abeta_const hM β₀]

/-- `thm:simple_thm1`, stacksvd half. With `θ_i = θ₀` and `c_i = c₀` the threshold is
`θ₀⁴ > c₀/M` and the limit is the paper's `1 - (c₀+θ₀²)/(Mθ₀⁴+θ₀²)`. -/
theorem simple_thm1_stacksvd (hM : 0 < M) {θ₀ c₀ : ℝ} (hc : 0 < c₀) :
    stackSVDLimit (fun _ : Fin M => θ₀) (fun _ : Fin M => c₀)
      = if c₀ < (M : ℝ) * θ₀ ^ 4 then
          1 - (c₀ + θ₀ ^ 2) / ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2) else 0 := by
  have hMpos : (0 : ℝ) < M := by exact_mod_cast hM
  have hT : ∑ _i : Fin M, θ₀ ^ 2 = (M : ℝ) * θ₀ ^ 2 := by
    rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  have hC : ∑ _i : Fin M, c₀ = (M : ℝ) * c₀ := by
    rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  simp only [stackSVDLimit, gt_iff_lt]
  rw [hT, hC]
  have hguard : (M : ℝ) * c₀ < ((M : ℝ) * θ₀ ^ 2) ^ 2 ↔ c₀ < (M : ℝ) * θ₀ ^ 4 := by
    constructor
    · intro h; nlinarith
    · intro h; nlinarith
  split_ifs with h1 h2 h2
  · have hθ4 : (0 : ℝ) < θ₀ ^ 4 := by nlinarith
    have hθ2 : (0 : ℝ) < θ₀ ^ 2 := by nlinarith [sq_nonneg θ₀]
    have hD : (0 : ℝ) < (M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2 := by nlinarith
    have hDne : ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2) ≠ 0 := hD.ne'
    have hT0 : (0 : ℝ) < (M : ℝ) * θ₀ ^ 2 := by nlinarith
    have hRHS : 1 - (c₀ + θ₀ ^ 2) / ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2)
        = ((M : ℝ) * θ₀ ^ 4 - c₀) / ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2) := by
      rw [show (M : ℝ) * θ₀ ^ 4 - c₀
        = ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2) - (c₀ + θ₀ ^ 2) by ring, sub_div, div_self hDne]
    rw [hRHS, div_eq_div_iff (mul_pos hT0 (by linarith)).ne' hDne]
    ring
  · exact absurd (hguard.mp h1) h2
  · exact absurd (hguard.mpr h2) h1
  · rfl

/-- `thm:simple_thm1`, svdstack half. With `θ_i = θ₀` and `c_i = c₀` the threshold is
`θ₀⁴ > c₀` and the limit is the paper's `1 - (c₀+θ₀²)/(Mθ₀⁴+θ₀²-(M-1)c₀)`. -/
theorem simple_thm1_svdstack (hM : 0 < M) {θ₀ c₀ : ℝ} (hc : 0 < c₀) :
    svdstackLimit (fun _ : Fin M => beta θ₀ c₀)
      = if c₀ < θ₀ ^ 4 then
          1 - (c₀ + θ₀ ^ 2) / ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2 - ((M : ℝ) - 1) * c₀) else 0 := by
  have hMpos : (0 : ℝ) < M := by exact_mod_cast hM
  rw [svdstackLimit_const hM (beta θ₀ c₀), beta_sq]
  simp only [betaSq, gt_iff_lt]
  split_ifs with h
  · have hθ4 : c₀ < θ₀ ^ 4 := h
    have hden : (0 : ℝ) < θ₀ ^ 4 + θ₀ ^ 2 := by nlinarith [sq_nonneg θ₀]
    have hdenne : (θ₀ ^ 4 + θ₀ ^ 2) ≠ 0 := hden.ne'
    have hD : (0 : ℝ) < (M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2 - ((M : ℝ) - 1) * c₀ := by
      nlinarith [mul_pos hMpos (sub_pos.mpr hθ4), sq_nonneg θ₀]
    have hDne : ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2 - ((M : ℝ) - 1) * c₀) ≠ 0 := hD.ne'
    have hsplit : 1 + ((M : ℝ) - 1) * ((θ₀ ^ 4 - c₀) / (θ₀ ^ 4 + θ₀ ^ 2))
        = ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2 - ((M : ℝ) - 1) * c₀) / (θ₀ ^ 4 + θ₀ ^ 2) := by
      rw [show (M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2 - ((M : ℝ) - 1) * c₀
        = (θ₀ ^ 4 + θ₀ ^ 2) + ((M : ℝ) - 1) * (θ₀ ^ 4 - c₀) by ring, add_div, div_self hdenne,
        mul_div_assoc]
    have hRHS : 1 - (c₀ + θ₀ ^ 2) / ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2 - ((M : ℝ) - 1) * c₀)
        = (M : ℝ) * (θ₀ ^ 4 - c₀) / ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2 - ((M : ℝ) - 1) * c₀) := by
      rw [show (M : ℝ) * (θ₀ ^ 4 - c₀)
        = ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2 - ((M : ℝ) - 1) * c₀) - (c₀ + θ₀ ^ 2) by ring, sub_div,
        div_self hDne]
    rw [hsplit, hRHS, div_eq_div_iff (div_pos hD hden).ne' hDne]
    field_simp
  · simp

end EqualTables

/-! ### The strict halves of `prop:dominance`

`main_paper.tex:633`: optimally weighted stacksvd "provid[es] strict improvement above its
recovery threshold when `θ_i²/c_i` is not constant across `i`, and when at least two `θ_i`
are nonzero (respectively)". The recovery threshold is `∑_i θ_i⁴/c_i > 1`; the paper leaves
the two other conditions in words, and they become `∃ i j, θ_i² c_j ≠ θ_j² c_i` and
`∃ i j, i ≠ j ∧ θ_i ≠ 0 ∧ θ_j ≠ 0`. Both are tight: with `θ_i²/c_i` constant the optimal
weights are the unit weights, and with one nonzero `θ_i` the two estimators agree.
-/

section StrictDominance

variable {M : ℕ}

/-- `x²/y` splits into its tangent line at `u` and a square: this is the equality case of the
Engel (Sedrakyan) form of Cauchy-Schwarz made explicit. -/
private theorem sq_div_eq_tangent_add {x y u : ℝ} (hy : y ≠ 0) :
    x ^ 2 / y = (2 * u * x - u ^ 2 * y) + (x - u * y) ^ 2 / y := by
  field_simp
  ring

/-- Strict Engel (Sedrakyan) form of Cauchy-Schwarz: `(∑ a)²/(∑ b) < ∑ a²/b` as soon as two
of the ratios `a_k/b_k` differ. Mathlib has the non-strict form
(`Finset.sq_sum_div_le_sum_sq_div`) but not this one. -/
theorem sq_sum_div_lt_sum_sq_div {ι : Type*} {s : Finset ι} {a b : ι → ℝ}
    (hb : ∀ k ∈ s, 0 < b k) {i j : ι} (hi : i ∈ s) (hj : j ∈ s)
    (hne : a i * b j ≠ a j * b i) :
    (∑ k ∈ s, a k) ^ 2 / (∑ k ∈ s, b k) < ∑ k ∈ s, a k ^ 2 / b k := by
  have hB : 0 < ∑ k ∈ s, b k := Finset.sum_pos hb ⟨i, hi⟩
  have hBne : (∑ k ∈ s, b k) ≠ 0 := hB.ne'
  set t : ℝ := (∑ k ∈ s, a k) / (∑ k ∈ s, b k) with htdef
  have hsplit : ∑ k ∈ s, a k ^ 2 / b k
      = ∑ k ∈ s, ((2 * t * a k - t ^ 2 * b k) + (a k - t * b k) ^ 2 / b k) :=
    Finset.sum_congr rfl fun k hk => sq_div_eq_tangent_add (hb k hk).ne'
  have hlin : ∑ k ∈ s, (2 * t * a k - t ^ 2 * b k)
      = (∑ k ∈ s, a k) ^ 2 / (∑ k ∈ s, b k) := by
    rw [Finset.sum_sub_distrib, ← Finset.mul_sum, ← Finset.mul_sum, htdef]
    field_simp
    ring
  have hpos : 0 < ∑ k ∈ s, (a k - t * b k) ^ 2 / b k := by
    have hnn : ∀ k ∈ s, 0 ≤ (a k - t * b k) ^ 2 / b k := fun k hk =>
      div_nonneg (sq_nonneg _) (hb k hk).le
    have hex : a i - t * b i ≠ 0 ∨ a j - t * b j ≠ 0 := by
      by_contra hcon
      simp only [not_or, not_not] at hcon
      exact hne (by rw [show a i = t * b i by linarith [hcon.1],
        show a j = t * b j by linarith [hcon.2]]; ring)
    rcases hex with hx | hx
    · exact Finset.sum_pos' hnn ⟨i, hi, div_pos (by positivity) (hb i hi)⟩
    · exact Finset.sum_pos' hnn ⟨j, hj, div_pos (by positivity) (hb j hj)⟩
  rw [hsplit, Finset.sum_add_distrib, hlin]
  linarith

/-- The strict companion of `le_stackSVDLimitW`. -/
theorem lt_stackSVDLimitW {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hthr : 1 < ∑ i, θ i ^ 4 / c i) {y : ℝ} (hy : y ∈ Set.Icc (0 : ℝ) 1)
    (hgy : 1 < gW θ c y) : y < stackSVDLimitW θ c := by
  obtain ⟨hx, hgx⟩ := stackSVDLimitW_spec (existsUnique_root hc hthr)
  by_contra hle
  rw [not_lt] at hle
  rcases eq_or_lt_of_le hle with heq | hlt
  · rw [← heq] at hgy
    linarith
  · have := gW_strictAntiOn hc (exists_ne_zero_of_thr hthr)
      (Set.Ioo_subset_Icc_self hx) hy hlt
    linarith

/-- The strict form of `dominance_term`. The slack of `eq:dominance_termwise_ineq` is
`c(θ²+c)(s - svdTerm θ c)`, so the inequality is strict exactly when the table is not the
only one carrying signal. -/
theorem dominance_term_lt {θ c s : ℝ} (hc : 0 < c) (hs : 0 < s) (hθ : θ ≠ 0)
    (hlt : svdTerm θ c < s) :
    svdTerm θ c / s < θ ^ 4 * (1 - s / (s + 1)) / (c + s / (s + 1) * θ ^ 2) := by
  have h1 : (0 : ℝ) < s + 1 := by linarith
  have hd : (0 : ℝ) < (s + 1) * c + s * θ ^ 2 :=
    add_pos_of_pos_of_nonneg (mul_pos h1 hc) (mul_nonneg hs.le (sq_nonneg _))
  have hθ2 : (0 : ℝ) < θ ^ 2 := by
    rcases lt_trichotomy θ 0 with h | h | h
    · nlinarith
    · exact absurd h hθ
    · nlinarith
  rw [gW_term_at_div hc hs.le]
  unfold svdTerm at hlt ⊢
  split_ifs with h
  · have hy : (0 : ℝ) < θ ^ 2 + c := by linarith
    rw [if_pos h] at hlt
    have hlt' : θ ^ 4 - c < s * (θ ^ 2 + c) := by
      rw [div_lt_iff₀ hy] at hlt
      linarith
    rw [div_div, div_lt_div_iff₀ (by positivity) hd]
    nlinarith [mul_pos hc (by linarith : (0 : ℝ) < s * (θ ^ 2 + c) - (θ ^ 4 - c))]
  · rw [zero_div]
    have h4 : (0 : ℝ) < θ ^ 4 := by nlinarith
    exact div_pos h4 hd

/-- `prop:dominance`, strict first half: above the recovery threshold, optimally weighted
stacksvd strictly beats optimally weighted svdstack as soon as two tables carry signal. With
one table the two estimators agree, so the hypothesis is tight. -/
theorem svdstackOpt_lt_stackSVDLimitW {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hthr : 1 < ∑ i, θ i ^ 4 / c i) (htwo : ∃ i j, i ≠ j ∧ θ i ≠ 0 ∧ θ j ≠ 0) :
    svdstackLimitOpt (fun i => beta (θ i) (c i)) < stackSVDLimitW θ c := by
  have hLpos : 0 < stackSVDLimitW θ c := (stackSVDLimitW_spec (existsUnique_root hc hthr)).1.1
  set s := Sval (fun i => beta (θ i) (c i)) with hsdef
  have hSval : s = ∑ i, svdTerm (θ i) (c i) := by rw [hsdef, Sval_beta θ c hc]
  have hs0 : 0 ≤ s := Sval_beta_nonneg θ c hc
  rcases eq_or_lt_of_le hs0 with hz | hpos
  · rw [svdstackLimitOpt, ← hsdef, ← hz]
    simpa using hLpos
  · obtain ⟨i1, j1, hij, hθi, hθj⟩ := htwo
    have hterms : ∀ k : Fin M, 0 ≤ svdTerm (θ k) (c k) := fun k => svdTerm_nonneg (hc k)
    have hidx : ∃ k : Fin M, θ k ≠ 0 ∧ svdTerm (θ k) (c k) < s := by
      by_cases hzj : svdTerm (θ j1) (c j1) = 0
      · exact ⟨j1, hθj, by rw [hzj]; exact hpos⟩
      · refine ⟨i1, hθi, ?_⟩
        have hjpos : 0 < svdTerm (θ j1) (c j1) :=
          lt_of_le_of_ne (hterms j1) (Ne.symm hzj)
        have hpair := Finset.sum_le_sum_of_subset_of_nonneg
          (Finset.subset_univ ({i1, j1} : Finset (Fin M)))
          (fun k _ _ => hterms k)
        rw [Finset.sum_pair hij] at hpair
        rw [hSval]
        linarith
    obtain ⟨k0, hk0θ, hk0lt⟩ := hidx
    have hmem : s / (s + 1) ∈ Set.Icc (0 : ℝ) 1 := by
      constructor
      · positivity
      · rw [div_le_one (by linarith)]
        linarith
    have hgoal : 1 < gW θ c (s / (s + 1)) := by
      have hterm : ∀ i ∈ (Finset.univ : Finset (Fin M)), svdTerm (θ i) (c i) / s
          ≤ θ i ^ 4 * (1 - s / (s + 1)) / (c i + s / (s + 1) * θ i ^ 2) := by
        intro i _
        refine dominance_term (hc i) hpos ?_
        rw [hSval]
        exact Finset.single_le_sum (fun j _ => svdTerm_nonneg (hc j)) (Finset.mem_univ i)
      have hstrict : svdTerm (θ k0) (c k0) / s
          < θ k0 ^ 4 * (1 - s / (s + 1)) / (c k0 + s / (s + 1) * θ k0 ^ 2) :=
        dominance_term_lt (hc k0) hpos hk0θ hk0lt
      calc (1 : ℝ) = (∑ i, svdTerm (θ i) (c i)) / s := by rw [← hSval]; field_simp
        _ = ∑ i, svdTerm (θ i) (c i) / s := Finset.sum_div _ _ _
        _ < gW θ c (s / (s + 1)) :=
            Finset.sum_lt_sum hterm ⟨k0, Finset.mem_univ _, hstrict⟩
    exact lt_stackSVDLimitW hc hthr hmem hgoal

/-- `prop:dominance`, strict second half: above the recovery threshold, optimally weighted
stacksvd strictly beats unweighted stacksvd as soon as `θ_i²/c_i` is not constant. With
`θ_i²/c_i` constant the unit weights are optimal, so the hypothesis is tight. -/
theorem stackSVDLimit_lt_stackSVDLimitW {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hthr : 1 < ∑ i, θ i ^ 4 / c i) (hnc : ∃ i j, θ i ^ 2 * c j ≠ θ j ^ 2 * c i) :
    stackSVDLimit θ c < stackSVDLimitW θ c := by
  have hLpos : 0 < stackSVDLimitW θ c := (stackSVDLimitW_spec (existsUnique_root hc hthr)).1.1
  unfold stackSVDLimit
  split_ifs with h
  · set T := ∑ i, θ i ^ 2 with hTdef
    set C := ∑ i, c i with hCdef
    have hM : 0 < M := by
      rcases Nat.eq_zero_or_pos M with rfl | h'
      · rw [hTdef, hCdef] at h; simp at h
      · exact h'
    have hne : (Finset.univ : Finset (Fin M)).Nonempty :=
      Finset.univ_nonempty_iff.mpr ⟨⟨0, hM⟩⟩
    have hC : 0 < C := Finset.sum_pos (fun i _ => hc i) hne
    have hT0 : 0 ≤ T := Finset.sum_nonneg fun i _ => sq_nonneg _
    have hT : 0 < T := by nlinarith
    have hTT : 0 < T ^ 2 + T := by positivity
    set u := (T ^ 2 - C) / (T ^ 2 + T) with hudef
    have hu : u ∈ Set.Icc (0 : ℝ) 1 := by
      constructor
      · rw [hudef]
        exact div_nonneg (by linarith) hTT.le
      · rw [hudef, div_le_one hTT]
        linarith
    have hbpos : ∀ i ∈ (Finset.univ : Finset (Fin M)),
        0 < c i * (T ^ 2 + T) + θ i ^ 2 * (T ^ 2 - C) := fun i _ =>
      add_pos_of_pos_of_nonneg (mul_pos (hc i) hTT)
        (mul_nonneg (sq_nonneg _) (by linarith))
    have hsumb : ∑ i, (c i * (T ^ 2 + T) + θ i ^ 2 * (T ^ 2 - C)) = T ^ 2 * (T + C) := by
      rw [Finset.sum_add_distrib, ← Finset.sum_mul, ← Finset.sum_mul, ← hTdef, ← hCdef]
      ring
    have hterm : ∀ i : Fin M, θ i ^ 4 * (1 - u) / (c i + u * θ i ^ 2)
        = (T + C) * (θ i ^ 2) ^ 2 / (c i * (T ^ 2 + T) + θ i ^ 2 * (T ^ 2 - C)) := by
      intro i
      have hden : 0 < c i + u * θ i ^ 2 :=
        add_pos_of_pos_of_nonneg (hc i) (mul_nonneg hu.1 (sq_nonneg _))
      have hbi : 0 < c i * (T ^ 2 + T) + θ i ^ 2 * (T ^ 2 - C) := hbpos i (Finset.mem_univ i)
      rw [div_eq_div_iff hden.ne' hbi.ne', hudef]
      field_simp
      ring
    -- the strict Cauchy-Schwarz step
    obtain ⟨i0, j0, hne0⟩ := hnc
    have hcross : (fun i => θ i ^ 2) i0 * (fun i => c i * (T ^ 2 + T) + θ i ^ 2 * (T ^ 2 - C)) j0
        ≠ (fun i => θ i ^ 2) j0 * (fun i => c i * (T ^ 2 + T) + θ i ^ 2 * (T ^ 2 - C)) i0 := by
      simp only
      intro heq
      refine hne0 ?_
      have h0 : (T ^ 2 + T) * (θ i0 ^ 2 * c j0 - θ j0 ^ 2 * c i0) = 0 := by
        linear_combination heq
      rcases mul_eq_zero.mp h0 with h1 | h1
      · exact absurd h1 hTT.ne'
      · linarith
    have hcs := sq_sum_div_lt_sum_sq_div (s := (Finset.univ : Finset (Fin M)))
      (a := fun i => θ i ^ 2) (b := fun i => c i * (T ^ 2 + T) + θ i ^ 2 * (T ^ 2 - C))
      hbpos (Finset.mem_univ i0) (Finset.mem_univ j0) hcross
    rw [hsumb, ← hTdef] at hcs
    have hTC : 0 < T + C := by linarith
    have hfrac : T ^ 2 / (T ^ 2 * (T + C)) = 1 / (T + C) := by
      rw [mul_comm]
      field_simp
    rw [hfrac] at hcs
    have hgu : 1 < gW θ c u := by
      have hgeq : gW θ c u
          = (T + C) * ∑ i, (θ i ^ 2) ^ 2 / (c i * (T ^ 2 + T) + θ i ^ 2 * (T ^ 2 - C)) := by
        unfold gW
        rw [Finset.mul_sum]
        exact Finset.sum_congr rfl fun i _ => by rw [hterm i, mul_div_assoc]
      rw [hgeq]
      calc (1 : ℝ) = (T + C) * (1 / (T + C)) := by field_simp
        _ < (T + C) * ∑ i, (θ i ^ 2) ^ 2 / (c i * (T ^ 2 + T) + θ i ^ 2 * (T ^ 2 - C)) :=
            mul_lt_mul_of_pos_left hcs hTC
    have hres := lt_stackSVDLimitW hc hthr hu hgu
    rw [hudef] at hres
    rw [show T * (T + 1) = T ^ 2 + T by ring]
    exact hres
  · exact hLpos

end StrictDominance

end Scalars
end StackedSVD

namespace StackedSVD.Scalars

variable {M : ℕ}

/-! ### The best binary weighting: the maximum over subsets (`main_paper.tex:442`)

`cor.2` fixes one subset, `{i : θ_i⁴ > c_i}`. The paper then optimizes over every subset and
writes the limit as a maximum of the same ratio. `binaryStackSVDLimitMax` is that maximum.
Three facts follow: the maximum is attained (`exists_subset_eq_max`); it is at least the
unweighted stacksvd limit (`stackSVDLimit_le_binaryStackSVDLimitMax`, the paper's "trivially
always perform at least as well as unweighted stacksvd"); and the guard commutes with the
maximum (`binaryStackSVDLimitMax_eq_sup_ratio`, the paper's two-case display). The comparison
with the optimally weighted limit is `binaryStackSVDLimitMax_le_stackSVDLimitW` in
`StackSVDWeighted.lean`, where `Lw_binary` and `L_le_opt` live.
-/

/-- The `2^M` subsets of `[M]` form a nonempty `Finset`. -/
theorem powerset_univ_nonempty : (Finset.univ : Finset (Finset (Fin M))).Nonempty :=
  ⟨∅, Finset.mem_univ _⟩

/-- The limit of optimally binary-weighted stacksvd (`main_paper.tex:442`): the maximum of
`binaryStackSVDLimit` over every subset of the `M` tables. -/
noncomputable def binaryStackSVDLimitMax (θ c : Fin M → ℝ) : ℝ :=
  (Finset.univ : Finset (Finset (Fin M))).sup' powerset_univ_nonempty
    fun S => binaryStackSVDLimit S θ c

/-- The paper's ratio at a subset, with the guard removed (`main_paper.tex:442`, first
branch). `binaryStackSVDLimit` is this value when the subset is above the stacksvd threshold,
and `0` otherwise. -/
noncomputable def binaryStackSVDRatio (θ c : Fin M → ℝ) (S : Finset (Fin M)) : ℝ :=
  ((∑ i ∈ S, θ i ^ 2) ^ 2 - ∑ i ∈ S, c i) / ((∑ i ∈ S, θ i ^ 2) ^ 2 + ∑ i ∈ S, θ i ^ 2)

/-- Every subset is at most the maximum. -/
theorem binaryStackSVDLimit_le_max (S : Finset (Fin M)) (θ c : Fin M → ℝ) :
    binaryStackSVDLimit S θ c ≤ binaryStackSVDLimitMax θ c := by
  unfold binaryStackSVDLimitMax
  exact Finset.le_sup' (fun S => binaryStackSVDLimit S θ c) (Finset.mem_univ S)

/-- The maximum is attained: some subset realizes it. -/
theorem exists_subset_eq_max (θ c : Fin M → ℝ) :
    ∃ S : Finset (Fin M), binaryStackSVDLimitMax θ c = binaryStackSVDLimit S θ c := by
  obtain ⟨S, -, hS⟩ :=
    Finset.exists_mem_eq_sup' (powerset_univ_nonempty (M := M))
      fun S => binaryStackSVDLimit S θ c
  exact ⟨S, hS⟩

/-- A bound that holds for every subset holds for the maximum. -/
theorem binaryStackSVDLimitMax_le {θ c : Fin M → ℝ} {a : ℝ}
    (h : ∀ S : Finset (Fin M), binaryStackSVDLimit S θ c ≤ a) :
    binaryStackSVDLimitMax θ c ≤ a := by
  unfold binaryStackSVDLimitMax
  exact Finset.sup'_le _ _ fun S _ => h S

@[simp] theorem binaryStackSVDLimit_empty (θ c : Fin M → ℝ) :
    binaryStackSVDLimit (∅ : Finset (Fin M)) θ c = 0 := by
  unfold binaryStackSVDLimit
  simp

@[simp] theorem binaryStackSVDRatio_empty (θ c : Fin M → ℝ) :
    binaryStackSVDRatio θ c (∅ : Finset (Fin M)) = 0 := by
  unfold binaryStackSVDRatio
  simp

/-- Discarding every table gives `0`, so the maximum is never negative. -/
theorem binaryStackSVDLimitMax_nonneg (θ c : Fin M → ℝ) : 0 ≤ binaryStackSVDLimitMax θ c := by
  have h := binaryStackSVDLimit_le_max (∅ : Finset (Fin M)) θ c
  rwa [binaryStackSVDLimit_empty] at h

/-- Every binary-weighted value is nonnegative: above the guard both the numerator and the
denominator are nonnegative, below it the value is `0`. -/
theorem binaryStackSVDLimit_nonneg (S : Finset (Fin M)) (θ c : Fin M → ℝ) :
    0 ≤ binaryStackSVDLimit S θ c := by
  unfold binaryStackSVDLimit
  by_cases h : (∑ i ∈ S, θ i ^ 2) ^ 2 > ∑ i ∈ S, c i
  · rw [if_pos h]
    have hT : (0 : ℝ) ≤ ∑ i ∈ S, θ i ^ 2 := Finset.sum_nonneg fun i _ => sq_nonneg _
    exact div_nonneg (by linarith) (by positivity)
  · rw [if_neg h]

/-- The maximum over subsets is attained at a **nonempty** subset. The empty subset gives `0`
and every value is nonnegative, so a singleton attains the maximum whenever the empty subset
does. This is what makes the model half of `cor.2` (which needs `NeZero S.card`) apply at the
maximizing subset. -/
theorem exists_nonempty_eq_max [NeZero M] (θ c : Fin M → ℝ) :
    ∃ S : Finset (Fin M), S.Nonempty ∧
      binaryStackSVDLimitMax θ c = binaryStackSVDLimit S θ c := by
  obtain ⟨S, hS⟩ := exists_subset_eq_max θ c
  rcases S.eq_empty_or_nonempty with rfl | hne
  · refine ⟨{(0 : Fin M)}, Finset.singleton_nonempty _, ?_⟩
    have h0 : binaryStackSVDLimitMax θ c = 0 := by rw [hS, binaryStackSVDLimit_empty]
    have hle := binaryStackSVDLimit_le_max ({(0 : Fin M)}) θ c
    have hnn := binaryStackSVDLimit_nonneg ({(0 : Fin M)}) θ c
    rw [h0] at hle ⊢
    linarith
  · exact ⟨S, hne, hS⟩

/-- Keeping every table is unweighted stacksvd. -/
theorem binaryStackSVDLimit_univ (θ c : Fin M → ℝ) :
    binaryStackSVDLimit Finset.univ θ c = stackSVDLimit θ c := by
  unfold binaryStackSVDLimit stackSVDLimit
  split_ifs with h
  · rw [show (∑ i, θ i ^ 2) * ((∑ i, θ i ^ 2) + 1) = (∑ i, θ i ^ 2) ^ 2 + ∑ i, θ i ^ 2 by ring]
  · rfl

/-- `main_paper.tex:442`: binary weighting "will trivially always perform at least as well as
unweighted stacksvd". -/
theorem stackSVDLimit_le_binaryStackSVDLimitMax (θ c : Fin M → ℝ) :
    stackSVDLimit θ c ≤ binaryStackSVDLimitMax θ c := by
  rw [← binaryStackSVDLimit_univ]
  exact binaryStackSVDLimit_le_max _ _ _

/-- Below the guard the raw ratio is nonpositive, so it never exceeds the guarded value.
`T_S = ∑_S θ_i² ≥ 0` makes the denominator nonnegative. -/
theorem binaryStackSVDRatio_le (θ c : Fin M → ℝ) (S : Finset (Fin M)) :
    binaryStackSVDRatio θ c S ≤ binaryStackSVDLimit S θ c := by
  have hT : 0 ≤ ∑ i ∈ S, θ i ^ 2 := Finset.sum_nonneg fun i _ => sq_nonneg _
  unfold binaryStackSVDRatio binaryStackSVDLimit
  split_ifs with h
  · exact le_rfl
  · exact div_nonpos_of_nonpos_of_nonneg (by linarith [not_lt.mp h]) (by positivity)

/-- `main_paper.tex:442`, the two-case display. The guard commutes with the maximum: with the
Lean convention `x / 0 = 0`, the maximum of the guarded values is the maximum of the paper's
raw ratios. The paper's second branch is the empty subset, which contributes `0` to both
sides. No hypothesis is needed: below the guard the ratio is `≤ 0`. -/
theorem binaryStackSVDLimitMax_eq_sup_ratio (θ c : Fin M → ℝ) :
    binaryStackSVDLimitMax θ c =
      (Finset.univ : Finset (Finset (Fin M))).sup' powerset_univ_nonempty
        (binaryStackSVDRatio θ c) := by
  refine le_antisymm (binaryStackSVDLimitMax_le fun S => ?_) ?_
  · unfold binaryStackSVDLimit
    split_ifs with h
    · exact Finset.le_sup' (binaryStackSVDRatio θ c) (Finset.mem_univ S)
    · calc (0:ℝ) = binaryStackSVDRatio θ c (∅ : Finset (Fin M)) :=
            (binaryStackSVDRatio_empty θ c).symm
        _ ≤ _ := Finset.le_sup' (binaryStackSVDRatio θ c) (Finset.mem_univ _)
  · exact Finset.sup'_le _ _ fun S _ =>
      (binaryStackSVDRatio_le θ c S).trans (binaryStackSVDLimit_le_max S θ c)

/-- `main_paper.tex:443`, "this improves not only the performance of stacksvd, but also the
detectability threshold": one subset above the stacksvd threshold makes the maximum positive,
whatever the full collection does. -/
theorem binaryStackSVDLimitMax_pos {θ c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    {S : Finset (Fin M)} (hne : S.Nonempty) (h : ∑ i ∈ S, c i < (∑ i ∈ S, θ i ^ 2) ^ 2) :
    0 < binaryStackSVDLimitMax θ c := by
  refine lt_of_lt_of_le ?_ (binaryStackSVDLimit_le_max S θ c)
  unfold binaryStackSVDLimit
  rw [if_pos h]
  have hC : 0 < ∑ i ∈ S, c i := Finset.sum_pos (fun i _ => hc i) hne
  have hT : 0 ≤ ∑ i ∈ S, θ i ^ 2 := Finset.sum_nonneg fun i _ => sq_nonneg _
  have hTpos : 0 < ∑ i ∈ S, θ i ^ 2 := by nlinarith
  exact div_pos (by linarith) (by nlinarith)

end StackedSVD.Scalars

namespace StackedSVD.Scalars

variable {M : ℕ}

/-! ### `prop:binarystacksvd_inadmissable` (`main_paper.tex:660`)

The paper builds one instance: `θ_i = 1` for every table and `c_j = 2j - 1` (one-based `j`),
which is `c i = 2i + 1` on the zero-based `Fin M`. On it,

1. every subset sits at or below the stacksvd threshold, since `|S|² ≤ ∑_{i ∈ S} (2i+1)`
   with equality exactly on the prefixes, so optimally binary-weighted stacksvd is `0`;
2. every table sits at its own threshold, `θ_i⁴ = 1 ≤ c_i`, so optimally weighted svdstack
   is `0`; unweighted stacksvd and unweighted svdstack are `0` for the same two reasons;
3. optimally weighted stacksvd is above `1 - ε` as soon as `log M + γ ≥ 2/ε`, because
   `gW` at `x = 1 - ε` is `(ε/2) ∑_{j=1}^M 1/(j - ε/2) ≥ (ε/2) H_M > (ε/2)(log M + γ)`.
-/

/-- The instance of `prop:binarystacksvd_inadmissable`: every `θ_i` is `1`. -/
def inadTheta (M : ℕ) : Fin M → ℝ := fun _ => 1

/-- The instance of `prop:binarystacksvd_inadmissable`: `c_j = 2j - 1` in the paper's
one-based index, which is `c i = 2i + 1` on the zero-based `Fin M`. -/
def inadC (M : ℕ) : Fin M → ℝ := fun i => 2 * (i.val : ℝ) + 1

theorem inadC_pos (M : ℕ) : ∀ i, 0 < inadC M i := by
  intro i
  have hi : (0:ℝ) ≤ (i.val : ℝ) := Nat.cast_nonneg _
  unfold inadC
  linarith

/-! #### Step 1: every subset is below the stacksvd threshold -/

/-- The sum of `2i + 1` over any `n` distinct naturals is at least `n²`. Induction on `n`:
the largest element `m` of `S` has `m + 1 ≥ #S`, so `(k+1)² = k² + (2k+1) ≤ ∑_{S \ m} + (2m+1)`.
Equality holds exactly for the prefixes `{0, ..., n-1}`. -/
private theorem card_sq_le_sum_odd :
    ∀ (n : ℕ) (S : Finset ℕ), S.card = n → n ^ 2 ≤ ∑ i ∈ S, (2 * i + 1) := by
  intro n
  induction n with
  | zero => intro S _; simp
  | succ k ih =>
    intro S hS
    have hne : S.Nonempty := by
      rw [← Finset.card_pos, hS]; omega
    have hmem : S.max' hne ∈ S := S.max'_mem hne
    have hsub : S ⊆ Finset.range (S.max' hne + 1) := fun x hx =>
      Finset.mem_range.mpr (Nat.lt_succ_of_le (S.le_max' x hx))
    have hcard : S.card ≤ S.max' hne + 1 := by
      simpa using Finset.card_le_card hsub
    have herase : (S.erase (S.max' hne)).card = k := by
      rw [Finset.card_erase_of_mem hmem, hS]
      omega
    have hrec := ih _ herase
    have hsum : (2 * S.max' hne + 1) + ∑ i ∈ S.erase (S.max' hne), (2 * i + 1)
        = ∑ i ∈ S, (2 * i + 1) := Finset.add_sum_erase S (fun i => 2 * i + 1) hmem
    have hexp : (k + 1) ^ 2 = k ^ 2 + (2 * k + 1) := by ring
    rw [← hsum, hexp]
    omega

private theorem card_sq_le_sum_odd_fin (S : Finset (Fin M)) :
    S.card ^ 2 ≤ ∑ i ∈ S, (2 * i.val + 1) := by
  classical
  have hinj : Set.InjOn Fin.val (S : Set (Fin M)) := fun a _ b _ h => Fin.val_injective h
  have himg : ∑ x ∈ S.image Fin.val, (2 * x + 1) = ∑ i ∈ S, (2 * i.val + 1) :=
    Finset.sum_image hinj
  have hcard : (S.image Fin.val).card = S.card :=
    Finset.card_image_of_injective S Fin.val_injective
  rw [← himg, ← hcard]
  exact card_sq_le_sum_odd _ _ rfl

/-- On the instance, `(∑_{i ∈ S} θ_i²)² = |S|² ≤ ∑_{i ∈ S} c_i` for every subset `S`. -/
theorem sq_card_le_sum_inadC (S : Finset (Fin M)) :
    ((S.card : ℝ)) ^ 2 ≤ ∑ i ∈ S, inadC M i := by
  have h := card_sq_le_sum_odd_fin S
  have hcast : ((∑ i ∈ S, (2 * i.val + 1) : ℕ) : ℝ) = ∑ i ∈ S, inadC M i := by
    simp only [inadC]
    push_cast
    rfl
  rw [← hcast]
  exact_mod_cast h

theorem inad_sum_theta_sq (S : Finset (Fin M)) :
    ∑ i ∈ S, inadTheta M i ^ 2 = (S.card : ℝ) := by
  simp [inadTheta]

/-- Every binary weighting of the instance is below the stacksvd threshold. -/
theorem inad_binaryStackSVDLimit_eq_zero (S : Finset (Fin M)) :
    binaryStackSVDLimit S (inadTheta M) (inadC M) = 0 := by
  have hg : ¬ ((∑ i ∈ S, inadTheta M i ^ 2) ^ 2 > ∑ i ∈ S, inadC M i) := by
    rw [inad_sum_theta_sq]
    exact not_lt.mpr (sq_card_le_sum_inadC S)
  unfold binaryStackSVDLimit
  exact if_neg hg

/-- Optimally binary-weighted stacksvd fails on the instance. -/
theorem inad_binaryStackSVDLimitMax_eq_zero (M : ℕ) :
    binaryStackSVDLimitMax (inadTheta M) (inadC M) = 0 :=
  le_antisymm
    (binaryStackSVDLimitMax_le fun S => le_of_eq (inad_binaryStackSVDLimit_eq_zero S))
    (binaryStackSVDLimitMax_nonneg _ _)

/-- Unweighted stacksvd fails on the instance. -/
theorem inad_stackSVDLimit_eq_zero (M : ℕ) :
    stackSVDLimit (inadTheta M) (inadC M) = 0 := by
  rw [← binaryStackSVDLimit_univ]
  exact inad_binaryStackSVDLimit_eq_zero _

/-! #### Step 2: every table is at its own detectability threshold -/

private theorem inad_le (M : ℕ) (i : Fin M) : inadTheta M i ^ 4 ≤ inadC M i := by
  have hi : (0:ℝ) ≤ (i.val : ℝ) := Nat.cast_nonneg _
  simp only [inadTheta, inadC, one_pow]
  linarith

private theorem inad_svdTerm_eq_zero (M : ℕ) (i : Fin M) :
    svdTerm (inadTheta M i) (inadC M i) = 0 := by
  unfold svdTerm
  exact if_neg (not_lt.mpr (inad_le M i))

private theorem inad_betaSq_eq_zero (M : ℕ) (i : Fin M) :
    betaSq (inadTheta M i) (inadC M i) = 0 := by
  unfold betaSq
  exact if_neg (not_lt.mpr (inad_le M i))

private theorem inad_beta_eq_zero (M : ℕ) (i : Fin M) :
    beta (inadTheta M i) (inadC M i) = 0 := by
  unfold beta
  rw [inad_betaSq_eq_zero, Real.sqrt_zero]

/-- Optimally weighted svdstack fails on the instance. -/
theorem inad_svdstackLimitOpt_eq_zero (M : ℕ) :
    svdstackLimitOpt (fun i => beta (inadTheta M i) (inadC M i)) = 0 := by
  have hS : Sval (fun i => beta (inadTheta M i) (inadC M i)) = 0 := by
    rw [Sval_beta _ _ (inadC_pos M)]
    exact Finset.sum_eq_zero fun i _ => inad_svdTerm_eq_zero M i
  unfold svdstackLimitOpt
  rw [hS]
  norm_num

/-- Unweighted svdstack fails on the instance. -/
theorem inad_svdstackLimit_eq_zero (M : ℕ) :
    svdstackLimit (fun i => beta (inadTheta M i) (inadC M i)) = 0 := by
  have hb : (fun i => beta (inadTheta M i) (inadC M i)) = (0 : Fin M → ℝ) :=
    funext fun i => inad_beta_eq_zero M i
  rw [hb]
  unfold svdstackLimit
  simp

/-! #### Step 3: optimally weighted stacksvd tends to 1 -/

private theorem sum_inv_eq_harmonic (M : ℕ) :
    ∑ i : Fin M, (1:ℝ) / ((i.val : ℝ) + 1) = (harmonic M : ℝ) := by
  rw [Fin.sum_univ_eq_sum_range (fun j : ℕ => (1:ℝ) / ((j : ℝ) + 1)) M]
  simp only [harmonic]
  push_cast
  exact Finset.sum_congr rfl fun j _ => by rw [one_div]

/-- `H_M > log M + γ` for `M ≥ 1`, the paper's `∑ 1/i ≥ ln M + γ`. Mathlib's
`eulerMascheroniSeq' M = H_M - log M` decreases to `γ`. -/
private theorem log_add_gamma_lt_harmonic {M : ℕ} (hM : 1 ≤ M) :
    Real.log M + Real.eulerMascheroniConstant < (harmonic M : ℝ) := by
  have h := Real.eulerMascheroniConstant_lt_eulerMascheroniSeq' M
  rw [Real.eulerMascheroniSeq', if_neg (show ¬ (M = 0) by omega)] at h
  linarith

/-- `gW` on the instance at `x = 1 - ε`: each summand is `ε/(2i + 2 - ε)`. -/
private theorem inad_gW_eq (M : ℕ) (ε : ℝ) :
    gW (inadTheta M) (inadC M) (1 - ε) = ∑ i : Fin M, ε / (2 * (i.val : ℝ) + 2 - ε) := by
  unfold gW
  refine Finset.sum_congr rfl fun i _ => ?_
  simp only [inadTheta, inadC]
  congr 1
  · ring
  · ring

/-- The paper's chain: `∑ ε/(2i+2-ε) ≥ (ε/2) H_M > (ε/2)(log M + γ) ≥ 1`. -/
private theorem one_lt_inad_gW {M : ℕ} (hM : 1 ≤ M) {ε : ℝ} (hε : 0 < ε) (hε1 : ε < 1)
    (hlog : 2 / ε ≤ Real.log M + Real.eulerMascheroniConstant) :
    1 < gW (inadTheta M) (inadC M) (1 - ε) := by
  have hlow : ∀ i : Fin M,
      ε / 2 * (1 / ((i.val : ℝ) + 1)) ≤ ε / (2 * (i.val : ℝ) + 2 - ε) := by
    intro i
    have hi : (0:ℝ) ≤ (i.val : ℝ) := Nat.cast_nonneg _
    have e1 : ε / 2 * (1 / ((i.val : ℝ) + 1)) = ε / (2 * (i.val : ℝ) + 2) := by
      rw [div_mul_div_comm]
      congr 1
      · ring
      · ring
    rw [e1]
    exact div_le_div_of_nonneg_left hε.le (by linarith) (by linarith)
  have hsum : ε / 2 * (harmonic M : ℝ) ≤ gW (inadTheta M) (inadC M) (1 - ε) := by
    rw [inad_gW_eq]
    calc ε / 2 * (harmonic M : ℝ)
        = ∑ i : Fin M, ε / 2 * (1 / ((i.val : ℝ) + 1)) := by
          rw [← sum_inv_eq_harmonic M, Finset.mul_sum]
      _ ≤ ∑ i : Fin M, ε / (2 * (i.val : ℝ) + 2 - ε) := Finset.sum_le_sum fun i _ => hlow i
  have hpos : (0:ℝ) < ε / 2 := by linarith
  have s1 : ε / 2 * (2 / ε) ≤ ε / 2 * (Real.log M + Real.eulerMascheroniConstant) :=
    mul_le_mul_of_nonneg_left hlog hpos.le
  have s2 : ε / 2 * (Real.log M + Real.eulerMascheroniConstant) < ε / 2 * (harmonic M : ℝ) :=
    mul_lt_mul_of_pos_left (log_add_gamma_lt_harmonic hM) hpos
  have s0 : ε / 2 * (2 / ε) = 1 := by field_simp
  linarith

/-- The instance is above the weighted stacksvd threshold. -/
private theorem inad_thr {M : ℕ} (hM : 1 ≤ M) {ε : ℝ} (hε : 0 < ε) (hε1 : ε < 1)
    (hg : 1 < gW (inadTheta M) (inadC M) (1 - ε)) :
    1 < ∑ i, inadTheta M i ^ 4 / inadC M i := by
  have hθ : ∃ i : Fin M, inadTheta M i ≠ 0 := ⟨⟨0, hM⟩, one_ne_zero⟩
  have hanti := gW_strictAntiOn (inadC_pos M) hθ
    (Set.left_mem_Icc.mpr zero_le_one) (Set.mem_Icc.mpr ⟨by linarith, by linarith⟩)
    (show (0:ℝ) < 1 - ε by linarith)
  rw [gW_zero] at hanti
  linarith

/-- Optimally weighted stacksvd on the instance, above `1 - ε` as soon as
`2/ε ≤ log M + γ`. -/
theorem inad_lt_stackSVDLimitW {M : ℕ} (hM : 1 ≤ M) {ε : ℝ} (hε : 0 < ε) (hε1 : ε < 1)
    (hlog : 2 / ε ≤ Real.log M + Real.eulerMascheroniConstant) :
    1 - ε < stackSVDLimitW (inadTheta M) (inadC M) := by
  have hg := one_lt_inad_gW hM hε hε1 hlog
  exact lt_stackSVDLimitW (inadC_pos M) (inad_thr hM hε hε1 hg)
    (Set.mem_Icc.mpr ⟨by linarith, by linarith⟩) hg

/-- `prop:binarystacksvd_inadmissable` (`main_paper.tex:660`): for every `ε ∈ (0,1)` and every
`M ≥ e^{-γ} exp(2/ε)`, the instance `θ_i = 1`, `c_i = 2i - 1` puts optimally weighted stacksvd
above `1 - ε` while optimally binary-weighted stacksvd and optimally weighted svdstack are
both `0`, and so are unweighted stacksvd and unweighted svdstack. -/
theorem binarystacksvd_inadmissable {ε : ℝ} (hε : 0 < ε) (hε1 : ε < 1) (M : ℕ)
    (hM : Real.exp (-Real.eulerMascheroniConstant) * Real.exp (2 / ε) ≤ M) :
    1 - ε < stackSVDLimitW (inadTheta M) (inadC M) ∧
      binaryStackSVDLimitMax (inadTheta M) (inadC M) = 0 ∧
      svdstackLimitOpt (fun i => beta (inadTheta M i) (inadC M i)) = 0 ∧
      stackSVDLimit (inadTheta M) (inadC M) = 0 ∧
      svdstackLimit (fun i => beta (inadTheta M i) (inadC M i)) = 0 := by
  have hexp : Real.exp (-Real.eulerMascheroniConstant) * Real.exp (2 / ε)
      = Real.exp (2 / ε - Real.eulerMascheroniConstant) := by
    rw [← Real.exp_add]
    congr 1
    ring
  rw [hexp] at hM
  have hγ : Real.eulerMascheroniConstant < 2 / 3 := Real.eulerMascheroniConstant_lt_two_thirds
  have h2ε : 2 < 2 / ε := by
    rw [lt_div_iff₀ hε]; linarith
  have hMlt : (1:ℝ) < M := by
    have h := Real.exp_lt_exp.mpr (show (0:ℝ) < 2 / ε - Real.eulerMascheroniConstant by linarith)
    rw [Real.exp_zero] at h
    linarith
  have hM1 : 1 ≤ M := by exact_mod_cast hMlt.le
  have hlog : 2 / ε ≤ Real.log M + Real.eulerMascheroniConstant := by
    have h0 : (0:ℝ) < M := by linarith
    have h := (Real.le_log_iff_exp_le h0).mpr hM
    linarith
  exact ⟨inad_lt_stackSVDLimitW hM1 hε hε1 hlog, inad_binaryStackSVDLimitMax_eq_zero M,
    inad_svdstackLimitOpt_eq_zero M, inad_stackSVDLimit_eq_zero M, inad_svdstackLimit_eq_zero M⟩

/-- `prop:binarystacksvd_inadmissable` at the paper's size `M = ⌈e^{-γ} exp(2/ε)⌉`. -/
theorem binarystacksvd_inadmissable_ceil {ε : ℝ} (hε : 0 < ε) (hε1 : ε < 1) :
    ∃ M : ℕ, M = ⌈Real.exp (-Real.eulerMascheroniConstant) * Real.exp (2 / ε)⌉₊ ∧
      1 - ε < stackSVDLimitW (inadTheta M) (inadC M) ∧
      binaryStackSVDLimitMax (inadTheta M) (inadC M) = 0 ∧
      svdstackLimitOpt (fun i => beta (inadTheta M i) (inadC M i)) = 0 ∧
      stackSVDLimit (inadTheta M) (inadC M) = 0 ∧
      svdstackLimit (fun i => beta (inadTheta M i) (inadC M i)) = 0 :=
  ⟨_, rfl, binarystacksvd_inadmissable hε hε1 _ (Nat.le_ceil _)⟩

end StackedSVD.Scalars
