/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.SingleWeight.Scalars
import StackedSVD.RankR.Example
import StackedSVD.RankR.General

/-!
# The suboptimality instance of `prop:singleweight_suboptimality`

Unit S5 of `notes/archive/singleweight_plan.md` section 4.3. The scalars only; the model and the
existence statement are `RankR/SingleWeight/Existence.lean`.

## The instance

`main_paper.tex:2171` to `:2194`. Two tables, `M = 2`, `r = 2`, one spike each
(`r_1 = r_2 = 1`), with

```
R_1 = (1, 0)ᵀ,   R_2 = (0, 1)ᵀ,   Θ_1 = Θ_2 = θ_0,   θ_0⁴ > c_0,   c_1 = c_2 = c_0.
```

Optimally weighted svdstack (which is unweighted svdstack here) reaches `2 β_0²` with
`β_0² = betaSq θ_0 c_0`. Every single weighting of stacksvd stays strictly below it. The
paper's closed form is the display at `main_paper.tex:2190` to `:2191`.

## The three branches

`secMat` is diagonal in this instance, so the two secular roots are the scalar roots
`γ_ℓ = w_ℓ²(1 + θ_0²)` with `z_ℓ = e_ℓ`. Each root is **detectable** when it lies above every
`w_i²` and its threshold sum of `assum:gen_rank_stacksvd_eig_sep` (`main_paper.tex:2106`) is
below `1`; that is `DetectableEx` below. `swLimitEx` is total in `w`: it is the sum of the two
terms when both roots are detectable, the one term when exactly one is, and `0` when neither
is. The paper covers the last two cases by the sentence at `main_paper.tex:2187`
("at least one component is undetectable and the performance is bounded above by `β_0²`") and
the sentence at `main_paper.tex:2194` (`w_1 = w_2`, where the two roots tie and the value is
`2 betaSq θ_0 (2 c_0)`; scan 2.4 of the plan measured the tied point to agree with the general
formula to 1.1e-16, so it is inside branch 1 and not a corner).

Numeric check (plan section 2.4): at `θ_0 = 1.6`, `c_0 = 1`, `w_1 = 1` the general formula and
this closed form agree to 2.2e-16 at `w_2² = 0.75, 0.90, 1.00`, with values `0.909678893`,
`0.988392985`, `0.999297753`, all below `2 β_0² = 1.218750`.

STATUS 2026-09-05: both theorems are proved.
-/

open scoped InnerProductSpace Matrix

namespace StackedSVD

namespace SingleWeight

/-! ### 1. The scalars of the instance -/

/-- The other index of `Fin 2`: `other 0 = 1` and `other 1 = 0`. It is the paper's `j ≠ ℓ`. -/
def other (l : Fin 2) : Fin 2 := if l = 0 then 1 else 0

/-- `γ_ℓ = w_ℓ²(1 + θ_0²)`, the `ℓ`-th secular root of the instance
(`main_paper.tex:2190`). -/
noncomputable def gammaEx (θ₀ : ℝ) (w : Fin 2 → ℝ) (l : Fin 2) : ℝ := w l ^ 2 * (1 + θ₀ ^ 2)

/-- The `ℓ`-th summand of the display at `main_paper.tex:2190` to `:2191`:

```
(θ_0⁴ - c_0 - c_0 w_j⁴ θ_0⁴/(w_ℓ²(1+θ_0²) - w_j²)²) / (θ_0²(1 + θ_0²)),      j = other ℓ.
```

It equals `β_0² - c_0 w_j⁴ θ_0⁴ / ((γ_ℓ - w_j²)² θ_0² (1 + θ_0²))`, because
`β_0² = (θ_0⁴ - c_0)/(θ_0²(1 + θ_0²))` shares the denominator. So each term is strictly below
`β_0²` whenever `w_j ≠ 0`, which is the route to `swLimitEx_lt` (plan scan 2.4, observation
1, checked to 3.6e-15 over 10 cases). -/
noncomputable def swTermEx (θ₀ c₀ : ℝ) (w : Fin 2 → ℝ) (l : Fin 2) : ℝ :=
  (θ₀ ^ 4 - c₀ - c₀ * w (other l) ^ 4 * θ₀ ^ 4 /
      (gammaEx θ₀ w l - w (other l) ^ 2) ^ 2) / (θ₀ ^ 2 * (1 + θ₀ ^ 2))

/-- Root `ℓ` of the instance is detectable: `γ_ℓ = w_ℓ²(1 + θ_0²)` lies above every `w_i²`,
and the threshold sum of `assum:gen_rank_stacksvd_eig_sep` (`main_paper.tex:2106`) is below
`1`. The two clauses are `IsSecularRoot`'s first clause and `EigSep.thresh` at `c = c_0`,
written on the instance. -/
def DetectableEx (θ₀ c₀ : ℝ) (w : Fin 2 → ℝ) (l : Fin 2) : Prop :=
  Scalars.wSqMax w < gammaEx θ₀ w l ∧
    ∑ i, c₀ * w i ^ 4 / (gammaEx θ₀ w l - w i ^ 2) ^ 2 < 1

open Classical in
/-- The closed form of the instance, total in `w` (`main_paper.tex:2190` to `:2191`). Three
branches, in this order.

1. Both roots detectable: the sum of the two terms of the paper's display.
2. Exactly one root detectable: that one term. The other component sits in the bulk and
   contributes nothing (`main_paper.tex:2187`).
3. Neither root detectable: `0`.

The `if` conditions are `DetectableEx θ₀ c₀ w 0` and `DetectableEx θ₀ c₀ w 1`, in that
order. -/
noncomputable def swLimitEx (θ₀ c₀ : ℝ) (w : Fin 2 → ℝ) : ℝ :=
  if DetectableEx θ₀ c₀ w 0 then
    if DetectableEx θ₀ c₀ w 1 then
      swTermEx θ₀ c₀ w 0 + swTermEx θ₀ c₀ w 1
    else
      swTermEx θ₀ c₀ w 0
  else
    if DetectableEx θ₀ c₀ w 1 then
      swTermEx θ₀ c₀ w 1
    else
      0

/-! ### 2. The strict inequality -/

/-- `w (other l)` is positive when both weights are. -/
private theorem other_pos (w : Fin 2 → ℝ) (h0 : 0 < w 0) (h1 : 0 < w 1) (l : Fin 2) :
    0 < w (other l) := by
  fin_cases l <;> simp only [other, Fin.zero_eta, Fin.mk_one, Fin.isValue] <;> assumption

/-- Each term of `swLimitEx` at `θ_0 = 8/5, c_0 = 1` is strictly below `39/64 = β_0²`, given
the other table's weight is positive: the penalty `w_j⁴θ_0⁴/(γ_ℓ - w_j²)²` is strictly
positive (plan scan 2.4, observation 1). -/
private theorem swTermEx_lt (w : Fin 2 → ℝ) (l : Fin 2) (hj : 0 < w (other l))
    (h : DetectableEx (8 / 5) 1 w l) : swTermEx (8 / 5) 1 w l < 39 / 64 := by
  have hsub : 0 < gammaEx (8 / 5) w l - w (other l) ^ 2 := by
    linarith [Scalars.le_wSqMax w (other l), h.1]
  have hpen : 0 < 1 * w (other l) ^ 4 * (8 / 5 : ℝ) ^ 4 /
      (gammaEx (8 / 5) w l - w (other l) ^ 2) ^ 2 := by positivity
  unfold swTermEx
  rw [div_lt_iff₀ (by norm_num : (0:ℝ) < (8 / 5 : ℝ) ^ 2 * (1 + (8 / 5 : ℝ) ^ 2))]
  nlinarith [hpen]

/-- The scalar half of **`prop:singleweight_suboptimality`** (`main_paper.tex:915`): at
`θ_0 = 8/5` and `c_0 = 1` every single weighting of stacksvd with both weights positive stays
strictly below `2 β_0² = 39/32`, the value of unweighted svdstack
(`main_paper.tex:2191`).

`8/5` keeps every constant rational: `θ_0² = 64/25`, `θ_0⁴ = 4096/625`, `β_0² = 39/64` and
`2 β_0² = 39/32`. The bound holds on all three branches of `swLimitEx`: on branch 1 each term
is `β_0²` minus a strictly positive penalty, on branch 2 the single term is below `β_0²`, and
on branch 3 the value is `0`. -/
theorem swLimitEx_lt : ∀ w : Fin 2 → ℝ, 0 < w 0 → 0 < w 1 →
    swLimitEx (8 / 5) 1 w < 2 * betaSq (8 / 5) 1 := by
  intro w h0 h1
  have hb : betaSq (8 / 5) 1 = 39 / 64 := by norm_num [betaSq]
  rw [hb]
  unfold swLimitEx
  by_cases hd0 : DetectableEx (8 / 5) 1 w 0 <;> by_cases hd1 : DetectableEx (8 / 5) 1 w 1
  · rw [if_pos hd0, if_pos hd1]
    linarith [swTermEx_lt w 0 (other_pos w h0 h1 0) hd0,
      swTermEx_lt w 1 (other_pos w h0 h1 1) hd1]
  · rw [if_pos hd0, if_neg hd1]
    linarith [swTermEx_lt w 0 (other_pos w h0 h1 0) hd0]
  · rw [if_neg hd0, if_pos hd1]
    linarith [swTermEx_lt w 1 (other_pos w h0 h1 1) hd1]
  · rw [if_neg hd0, if_neg hd1]
    norm_num

/-! ### 3. The tie to the general layer -/

/-- `sigMat` of the instance is the rank-one outer product supported on the single index `i`:
`R_i = e_i`, so `R_i Θ_i² R_iᵀ` is `θ_0²` at `(i, i)` and `0` elsewhere. -/
theorem sigMat_ex (θ₀ : ℝ) (i : Fin 2) :
    sigMat (rk := fun _ : Fin 2 => 1) (fun _ _ => θ₀) (Rone (RankR.Example.Rex 0)) i
      = Matrix.of fun k k' => if k = i ∧ k' = i then θ₀ ^ 2 else 0 := by
  ext k k'
  fin_cases i <;> fin_cases k <;> fin_cases k' <;>
    simp [sigMat, Rone, RankR.Example.Rex, Matrix.mul_apply]

/-- The quadratic form of `secDerivMat` at the standard basis vector `e_l` picks out the
`(l, l)` entry, a single term of the sum since `sigMat_ex` is supported on `i = l`. -/
theorem quad_ex (θ₀ : ℝ) (w : Fin 2 → ℝ) (γ : ℝ) (l : Fin 2) :
    WithLp.ofLp (EuclideanSpace.single l (1 : ℝ)) ⬝ᵥ
      (secDerivMat (rk := fun _ : Fin 2 => 1) (fun _ _ => θ₀) (Rone (RankR.Example.Rex 0)) w γ
        *ᵥ WithLp.ofLp (EuclideanSpace.single l (1 : ℝ)))
      = w l ^ 2 * θ₀ ^ 2 / (γ - w l ^ 2) ^ 2 := by
  simp only [secDerivMat, sigMat_ex]
  fin_cases l <;>
    simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two, Matrix.sum_apply] <;> ring

/-- `swTerm` of the general layer, read at the instance's root `γ_ℓ` and eigenvector `e_ℓ`, is
`swTermEx`: the shared denominator `θ_0²(1 + θ_0²)` cancels between the two once `γ_ℓ - w_ℓ²`
is rewritten as `w_ℓ² θ_0²` (module docstring, item `swTermEx`). -/
theorem swTerm_ex (θ₀ c₀ : ℝ) (w : Fin 2 → ℝ) (hθ : 0 < θ₀) (hw : ∀ i, 0 < w i)
    (l : Fin 2) (h : DetectableEx θ₀ c₀ w l) :
    swTerm (rk := fun _ : Fin 2 => 1) (fun _ _ => θ₀) (Rone (RankR.Example.Rex 0)) w
      (fun _ => c₀) (gammaEx θ₀ w l) (EuclideanSpace.single l 1) = swTermEx θ₀ c₀ w l := by
  have hsub0 : 0 < gammaEx θ₀ w l - w 0 ^ 2 := by linarith [Scalars.le_wSqMax w 0, h.1]
  have hsub1 : 0 < gammaEx θ₀ w l - w 1 ^ 2 := by linarith [Scalars.le_wSqMax w 1, h.1]
  unfold swTerm
  rw [quad_ex]
  fin_cases l <;>
    simp [Fin.sum_univ_two, other, swTermEx, gammaEx] at hsub0 hsub1 ⊢ <;>
    field_simp [(hw 0).ne', (hw 1).ne', hθ.ne', hsub0.ne', hsub1.ne'] <;>
    ring

/-- On the two-root branch the closed form of the instance is the general limit `swLimit` of
`main_paper.tex:2119`, read at `γ_ℓ = w_ℓ²(1 + θ_0²)`, `z_ℓ = e_ℓ`,
`R = Rone (RankR.Example.Rex 0)` (that is `R_1 = (1,0)ᵀ` and `R_2 = (0,1)ᵀ`, since
`sin 0 = 0` and `cos 0 = 1`), `Θ_i = θ_0` and `c_i = c_0`.

This is the step that lets `prop_gen_rank_stacksvd_singleweight` deliver the instance: it says
the general layer, specialized, is the paper's display at `main_paper.tex:2190` to `:2191`. -/
theorem swLimitEx_eq_swLimit (θ₀ c₀ : ℝ) (w : Fin 2 → ℝ) (hθ : 0 < θ₀) (hw : ∀ i, 0 < w i)
    (h0 : DetectableEx θ₀ c₀ w 0) (h1 : DetectableEx θ₀ c₀ w 1) :
    swLimitEx θ₀ c₀ w
      = swLimit (r := 2) (rk := fun _ : Fin 2 => 1) (fun _ _ => θ₀)
          (Rone (RankR.Example.Rex 0)) w (fun _ => c₀) (gammaEx θ₀ w)
          (fun l => EuclideanSpace.single l 1) := by
  unfold swLimitEx
  rw [if_pos h0, if_pos h1]
  unfold swLimit
  simp only [Fin.sum_univ_two]
  rw [swTerm_ex θ₀ c₀ w hθ hw 0 h0, swTerm_ex θ₀ c₀ w hθ hw 1 h1]

end SingleWeight

end StackedSVD
