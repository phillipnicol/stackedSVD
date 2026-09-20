/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RankR.RMT.Stack
import StackedSVD.RMT.TailShift

/-!
# The tail shift on `RankRStack`

Task A3a of the D32 campaign. The rank-`r` stack interface `RankRStack`
(`RankR/RMT/Stack.lean`) needs the same device `SpikedModel.shift` gives in
`RMT/TailShift.lean`: reindex `N ↦ N + k` so that a side condition of the shape
`∀ N, ns N = r + pp N`, `∀ N, 0 < pp N` can be produced from `Tendsto ns atTop atTop` alone,
at a large enough `k`. This closes the gap between a general regime hypothesis
(`Tendsto (fun N => (ns N : ℝ) / d N) atTop (𝓝 c)`) and the block-nonempty hypothesis every
downstream `RankRStack` theorem takes, so that task A3b can state the unconditional
`align_of_gaussian` once its input lands.

## Content

1. `RankRStack.shift`, `shift_GaussianNoise`, and the `rfl` bridges for the derived objects
   (`gram`, `core`, `coreEig`, `spikeVec`, `signalPart`).
2. `exists_shift_lt`: `Tendsto ns atTop atTop` gives a `k` past which `ns` clears any target
   rank `r`, in the `pp` form the chain of `RankRStack` theorems takes.
3. `tendsto_ns_atTop`: the regime limit `(ns N : ℝ) / d N → c > 0` with `d N → ∞` gives
   `ns N → ∞`, so `exists_shift_lt` applies without a separate hypothesis.
4. `shift_hdtop`, `shift_hns`: the regime hypotheses shift the same way `shift_Regime` does in
   `RMT/TailShift.lean`, and `tendstoInProb_normSq_specProjTop_of_shift`, the transfer of a
   limit of the projector overlap from the shifted stack to the stack.
5. A plumbing `example`: the universe-quantified conditional statement (`RankRStack` at an
   arbitrary probability space, with the `pp` side condition) gives the unconditional one.

STATUS: infrastructure only, not a paper theorem. See
`notes/archive/agent_reports/d32_A3a_shiftR.md`.
-/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix InnerProductSpace

universe u

namespace StackedSVD

namespace RankRStack

variable {Ω : ℕ → Type u} [∀ N, MeasurableSpace (Ω N)] {μ : ∀ N, Measure (Ω N)}
  {ns d : ℕ → ℕ} {r : ℕ}

/-! ### 1. The shifted stack -/

/-- The stack reindexed by `N ↦ N + k`: same `A`, `V`, `core`, evaluated at `N + k`, matching
`SpikedModel.shift` of `RMT/TailShift.lean`. -/
def shift (s : RankRStack μ ns d r) (k : ℕ) :
    RankRStack (fun N => μ (N + k)) (fun N => ns (N + k)) (fun N => d (N + k)) r where
  V := fun N => s.V (N + k)
  A := fun N => s.A (N + k)
  core := s.core
  Zu := fun N => s.Zu (N + k)
  E := fun N => s.E (N + k)
  X := fun N => s.X (N + k)
  hV := fun N => s.hV (N + k)
  hcore := fun N => s.hcore (N + k)
  hE := fun N => s.hE (N + k)
  hX := fun N => s.hX (N + k)
  hd := fun N => s.hd (N + k)
  hZmeas := fun N => s.hZmeas (N + k)

/-- The Gaussian noise law is stated index by index, so it shifts by evaluation. -/
theorem shift_GaussianNoise (s : RankRStack μ ns d r) (k : ℕ) (hG : s.GaussianNoise) :
    (s.shift k).GaussianNoise := fun N => hG (N + k)

section Bridges

variable (s : RankRStack μ ns d r) (k : ℕ)

@[simp] theorem shift_gram (N : ℕ) (ω : Ω (N + k)) :
    (s.shift k).gram N ω = s.gram (N + k) ω := rfl

@[simp] theorem shift_core : (s.shift k).core = s.core := rfl

@[simp] theorem shift_coreEig (j : Fin r) : (s.shift k).coreEig j = s.coreEig j := rfl

@[simp] theorem shift_spikeVec (j : Fin r) (N : ℕ) :
    (s.shift k).spikeVec j N = s.spikeVec j (N + k) := rfl

@[simp] theorem shift_signalPart (N : ℕ) :
    (s.shift k).signalPart N = s.signalPart (N + k) := rfl

end Bridges

/-! ### 2. The shift that clears a target rank -/

/-- `Tendsto ns atTop atTop` gives a shift after which `r < ns (N + k)` for every `N`, in the
`pp` form the theorems of this chain take. -/
theorem exists_shift_lt (r : ℕ) (hns : Tendsto ns atTop atTop) :
    ∃ k, ∃ pp : ℕ → ℕ, (∀ N, ns (N + k) = r + pp N) ∧ ∀ N, 0 < pp N := by
  obtain ⟨k, hk⟩ := eventually_atTop.mp (hns.eventually_ge_atTop (r + 1))
  refine ⟨k, fun N => ns (N + k) - r, fun N => ?_, fun N => ?_⟩ <;>
    · have h : r + 1 ≤ ns (N + k) := hk (N + k) (Nat.le_add_left k N)
      dsimp only
      omega

/-- `Tendsto (fun N => (ns N : ℝ) / d N) atTop (𝓝 c)` with `0 < c` and `Tendsto d atTop atTop`
give `Tendsto ns atTop atTop`: the regime limit alone produces the block-nonempty hypothesis
that `exists_shift_lt` needs. -/
theorem tendsto_ns_atTop {c : ℝ} (hc : 0 < c) (hdtop : Tendsto d atTop atTop)
    (hns : Tendsto (fun N => (ns N : ℝ) / d N) atTop (𝓝 c)) : Tendsto ns atTop atTop := by
  have hev : ∀ᶠ N in atTop, c / 2 < (ns N : ℝ) / (d N : ℝ) :=
    hns.eventually (eventually_gt_nhds (by linarith : c / 2 < c))
  have hdpos : ∀ᶠ N in atTop, (0 : ℕ) < d N := hdtop.eventually_gt_atTop 0
  have hdtopR : Tendsto (fun N => (d N : ℝ)) atTop atTop := tendsto_natCast_atTop_iff.mpr hdtop
  have hmul : Tendsto (fun N => c / 2 * (d N : ℝ)) atTop atTop :=
    hdtopR.const_mul_atTop (by linarith : (0 : ℝ) < c / 2)
  have hle : (fun N => c / 2 * (d N : ℝ)) ≤ᶠ[atTop] (fun N => (ns N : ℝ)) := by
    filter_upwards [hev, hdpos] with N hN hdN
    exact ((lt_div_iff₀ (by exact_mod_cast hdN)).mp hN).le
  exact tendsto_natCast_atTop_iff.mp (tendsto_atTop_mono' atTop hle hmul)

/-! ### 3. The shifted regime hypotheses -/

theorem shift_hdtop (k : ℕ) (hdtop : Tendsto d atTop atTop) :
    Tendsto (fun N => d (N + k)) atTop atTop :=
  (tendsto_add_atTop_iff_nat k).mpr hdtop

theorem shift_hns (k : ℕ) {c : ℝ}
    (hns : Tendsto (fun N => (ns N : ℝ) / d N) atTop (𝓝 c)) :
    Tendsto (fun N => (ns (N + k) : ℝ) / d (N + k)) atTop (𝓝 c) :=
  (tendsto_add_atTop_iff_nat (f := fun N => (ns N : ℝ) / (d N : ℝ)) k).mpr hns

/-! ### 4. The transfer -/

/-- **The transfer**: a limit in probability of the projector overlap for the shifted stack
gives the limit for the stack. The `rfl` lemmas of section 1 make the two functions agree, so
the shifted hypothesis is handed straight to `SpikedModel.tendstoInProb_of_shift`. -/
theorem tendstoInProb_normSq_specProjTop_of_shift (s : RankRStack μ ns d r) (k : ℕ)
    (j : Fin r) {a : ℝ}
    (h : TendstoInProb (fun N => μ (N + k)) (fun N ω =>
      ‖specProjTop ((s.shift k).gram N ω) ((s.shift k).isHermitian_gram N ω) r
        ((s.shift k).spikeVec j N)‖ ^ 2) a) :
    TendstoInProb μ (fun N ω =>
      ‖specProjTop (s.gram N ω) (s.isHermitian_gram N ω) r (s.spikeVec j N)‖ ^ 2) a :=
  SpikedModel.tendstoInProb_of_shift k h

/-! ### 5. Plumbing test

A conditional statement of the chain, given in its universe-quantified form (any probability
space, with the side condition `pp`), gives the unconditional statement (`RankRStack μ ns d r`
at the ambient `Ω`, under the regime limit alone). This is the shape task A3b's
`align_of_gaussian` will use: state the conditional theorem for an abstract `RankRStack` with
`hn`/`hp`, then discharge `hn`/`hp` by `exists_shift_lt` and transfer with
`tendstoInProb_normSq_specProjTop_of_shift`. -/

example [∀ N, IsProbabilityMeasure (μ N)]
    (H : ∀ {Ω' : ℕ → Type u} [∀ N, MeasurableSpace (Ω' N)] {μ' : ∀ N, Measure (Ω' N)}
      [∀ N, IsProbabilityMeasure (μ' N)] {ns' d' : ℕ → ℕ}
      (s : RankRStack μ' ns' d' r) (_hG : s.GaussianNoise) {c : ℝ} (_hc : 0 < c)
      (_hdtop : Tendsto d' atTop atTop)
      (_hns : Tendsto (fun N => (ns' N : ℝ) / d' N) atTop (𝓝 c))
      {pp : ℕ → ℕ} (_hn : ∀ N, ns' N = r + pp N) (_hp : ∀ N, 0 < pp N) (j : Fin r),
      TendstoInProb μ' (fun N ω => ‖specProjTop (s.gram N ω) (s.isHermitian_gram N ω) r
        (s.spikeVec j N)‖ ^ 2) (betaSq (Real.sqrt (s.coreEig j)) c))
    (s : RankRStack μ ns d r) (hG : s.GaussianNoise) {c : ℝ} (hc : 0 < c)
    (hdtop : Tendsto d atTop atTop)
    (hns : Tendsto (fun N => (ns N : ℝ) / d N) atTop (𝓝 c)) (j : Fin r) :
    TendstoInProb μ (fun N ω => ‖specProjTop (s.gram N ω) (s.isHermitian_gram N ω) r
      (s.spikeVec j N)‖ ^ 2) (betaSq (Real.sqrt (s.coreEig j)) c) := by
  obtain ⟨k, pp, hn, hp⟩ := exists_shift_lt r (tendsto_ns_atTop hc hdtop hns)
  have : ∀ N, IsProbabilityMeasure (μ (N + k)) := SpikedModel.isProbabilityMeasure_shift k
  exact s.tendstoInProb_normSq_specProjTop_of_shift k j
    (H (s.shift k) (s.shift_GaussianNoise k hG) hc (shift_hdtop k hdtop) (shift_hns k hns)
      hn hp j)

end RankRStack

end StackedSVD
