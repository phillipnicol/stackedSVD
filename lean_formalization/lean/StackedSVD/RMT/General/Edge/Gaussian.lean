/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.RMT.General.Edge.Arith
import StackedSVD.RMT.General.Edge.Trace
import StackedSVD.RMT.R3

/-! # Stage 3, unit G: the Gaussian trace bound

`∫ trace ((Y_GᵀY_G)^k) ≤ 2 d (d (bulkEdge c + ε))^k` at `k = k_N = ⌈C log d_N⌉`, for `N`
large (`notes/stage3_edge.md`, route C; `notes/STAGE3_CAMPAIGN.md`). The route is
`trace (W^k) ≤ d lamMax(W)^k`, `lamMax(W_G) ≤ ‖Z_G‖²/d`, `∫‖Z_G‖ ≤ √n + √d`
(`R3.integral_opNorm_le`) and the one-sided concentration `R3.measure_opNorm_ge_le`, through
one moment integral of the Gaussian tail. The statement is at the moment order of the
assembly only: at `k` of order `d` the tail term `(2μ)^(2k) exp(-κ² d/2)` of this route is
not small, and the bound itself fails there for small `ε` (Fable audit, 2026-09-10).

The `lamMax` step is done with the `l2` operator norm rather than with eigenvalues:
`trace (M) ≤ d ‖M‖` entry by entry, `‖(YᵀY)^k‖ ≤ ‖YᵀY‖^k` by submultiplicativity, and
`‖YᵀY‖ = ‖Y‖²` (`Matrix.l2_opNorm_conjTranspose_mul_self`). The tail of `‖Y‖^(2k)` above
`a = ∫‖Y‖ + κ√d` is paid for by one arithmetic-geometric mean split
`x^(2k) ≤ (θ/2) x^(4k) + 1/(2θ)` at `θ = exp(-κ² d/8)`, which needs only the fourth moment
`∫ ‖Y‖^(4k) ≤ 4 (nd)^(2k) (4k)!` of the Frobenius norm, not a layer cake. -/

open MeasureTheory ProbabilityTheory Filter Topology
open scoped Matrix Matrix.Norms.L2Operator ENNReal NNReal

namespace StackedSVD
namespace Edge

/-! ### One Gaussian entry: every moment, from the moment generating function -/

/-- `|x| ^ m ≤ m ! (e^x + e^{-x})`, the crude bound that gives every moment of one Gaussian. -/
private theorem abs_pow_le_factorial_mul (m : ℕ) (x : ℝ) :
    |x| ^ m ≤ (Nat.factorial m : ℝ) * (Real.exp x + Real.exp (-x)) := by
  have hfac : (0 : ℝ) < (Nat.factorial m : ℝ) := by exact_mod_cast m.factorial_pos
  have h1 : |x| ^ m / (Nat.factorial m : ℝ) ≤ Real.exp |x| :=
    Real.pow_div_factorial_le_exp _ (abs_nonneg x) m
  have h2 : Real.exp |x| ≤ Real.exp x + Real.exp (-x) := by
    rcases abs_cases x with ⟨hx, _⟩ | ⟨hx, _⟩
    · rw [hx]; have := Real.exp_pos (-x); linarith
    · rw [hx]; have := Real.exp_pos x; linarith
  rw [div_le_iff₀ hfac] at h1
  nlinarith [h1, h2, hfac]

/-- Every power of one standard Gaussian entry is integrable. -/
private theorem integrable_pow_gaussianReal (m : ℕ) :
    Integrable (fun x : ℝ => x ^ m) (gaussianReal 0 1) := by
  have he1 : Integrable (fun x : ℝ => Real.exp x) (gaussianReal 0 1) := by
    simpa using integrable_exp_mul_gaussianReal (μ := 0) (v := 1) 1
  have he2 : Integrable (fun x : ℝ => Real.exp (-x)) (gaussianReal 0 1) := by
    simpa using integrable_exp_mul_gaussianReal (μ := 0) (v := 1) (-1)
  refine Integrable.mono' ((he1.add he2).const_mul (Nat.factorial m : ℝ)) (by fun_prop) ?_
  filter_upwards with x
  rw [Real.norm_eq_abs, abs_pow]
  exact abs_pow_le_factorial_mul m x

/-- `∫ e^x ∂N(0,1) = e^{1/2}`, from the moment generating function. -/
private theorem integral_exp_gaussianReal :
    ∫ x : ℝ, Real.exp x ∂(gaussianReal 0 1) = Real.exp (1 / 2) := by
  have h : mgf id (gaussianReal 0 1) 1 = Real.exp (0 * 1 + (1 : ℝ≥0) * 1 ^ 2 / 2) := by
    rw [mgf_id_gaussianReal]
  simp only [mgf, id] at h
  simpa using h

/-- `∫ e^{-x} ∂N(0,1) = e^{1/2}`. -/
private theorem integral_exp_neg_gaussianReal :
    ∫ x : ℝ, Real.exp (-x) ∂(gaussianReal 0 1) = Real.exp (1 / 2) := by
  have h : mgf id (gaussianReal 0 1) (-1) = Real.exp (0 * (-1) + (1 : ℝ≥0) * (-1) ^ 2 / 2) := by
    rw [mgf_id_gaussianReal]
  simp only [mgf, id] at h
  simpa using h

/-- `∫ x^m ∂N(0,1) ≤ 4 m !`. The loss against the true `(m-1)!!` is a factor of order `m!`,
which the tail factor `exp(-κ² d / 8)` of the split below absorbs. -/
private theorem integral_pow_gaussianReal_le (m : ℕ) :
    ∫ x : ℝ, x ^ m ∂(gaussianReal 0 1) ≤ 4 * (Nat.factorial m : ℝ) := by
  have he1 : Integrable (fun x : ℝ => Real.exp x) (gaussianReal 0 1) := by
    simpa using integrable_exp_mul_gaussianReal (μ := 0) (v := 1) 1
  have he2 : Integrable (fun x : ℝ => Real.exp (-x)) (gaussianReal 0 1) := by
    simpa using integrable_exp_mul_gaussianReal (μ := 0) (v := 1) (-1)
  have hmono : ∫ x : ℝ, x ^ m ∂(gaussianReal 0 1)
      ≤ ∫ x : ℝ, (Nat.factorial m : ℝ) * (Real.exp x + Real.exp (-x)) ∂(gaussianReal 0 1) := by
    refine integral_mono (integrable_pow_gaussianReal m) ((he1.add he2).const_mul _) fun x => ?_
    exact le_trans (le_abs_self _) (by rw [abs_pow]; exact abs_pow_le_factorial_mul m x)
  rw [integral_const_mul, integral_add he1 he2, integral_exp_gaussianReal,
    integral_exp_neg_gaussianReal] at hmono
  have hexp : Real.exp (1 / 2) ≤ 2 := by
    have h1 : Real.exp (1 / 2) ^ 2 = Real.exp 1 := by
      rw [← Real.exp_nat_mul]; norm_num
    nlinarith [Real.exp_pos (1 / 2 : ℝ), Real.exp_one_lt_d9]
  have hfac : (0 : ℝ) ≤ (Nat.factorial m : ℝ) := by positivity
  nlinarith [hmono, hexp, hfac]

/-! ### One entry of a Gaussian matrix -/

/-- The matrix law is the nested product law, so an entry is a coordinate of a coordinate. -/
private theorem integral_entry_pow (n d m : ℕ) (i : Fin n) (j : Fin d) :
    ∫ Y : Matrix (Fin n) (Fin d) ℝ, (Y i j) ^ m ∂(gaussianMatrix n d)
      = ∫ x : ℝ, x ^ m ∂(gaussianReal 0 1) := by
  have hf : AEStronglyMeasurable (fun x : ℝ => x ^ m) (gaussianReal 0 1) := by fun_prop
  have e2 : ∫ row : Fin d → ℝ, (row j) ^ m ∂(Measure.pi fun _ : Fin d => gaussianReal 0 1)
      = ∫ x : ℝ, x ^ m ∂(gaussianReal 0 1) := PiLaw.integral_comp_coord _ hf j
  have hg : AEStronglyMeasurable (fun row : Fin d → ℝ => (row j) ^ m)
      (Measure.pi fun _ : Fin d => gaussianReal 0 1) :=
    ((measurable_pi_apply j).pow_const m).aestronglyMeasurable
  have e1 : ∫ Y : (Fin n) → (Fin d) → ℝ, (Y i j) ^ m
        ∂(Measure.pi fun _ : Fin n => Measure.pi fun _ : Fin d => gaussianReal 0 1)
      = ∫ row : Fin d → ℝ, (row j) ^ m ∂(Measure.pi fun _ : Fin d => gaussianReal 0 1) :=
    PiLaw.integral_comp_coord _ hg i
  have e0 : ∫ Y : Matrix (Fin n) (Fin d) ℝ, (Y i j) ^ m ∂(gaussianMatrix n d)
      = ∫ Y : (Fin n) → (Fin d) → ℝ, (Y i j) ^ m
          ∂(Measure.pi fun _ : Fin n => Measure.pi fun _ : Fin d => gaussianReal 0 1) := rfl
  rw [e0, e1, e2]

/-- Every power of one entry is integrable under the Gaussian matrix law. -/
private theorem integrable_entry_pow (n d m : ℕ) (i : Fin n) (j : Fin d) :
    Integrable (fun Y : Matrix (Fin n) (Fin d) ℝ => (Y i j) ^ m) (gaussianMatrix n d) := by
  have h1 : Integrable (fun row : Fin d → ℝ => (row j) ^ m)
      (Measure.pi fun _ : Fin d => gaussianReal 0 1) :=
    PiLaw.integrable_comp_coord _ (integrable_pow_gaussianReal m) j
  have h2 : Integrable (fun Y : (Fin n) → (Fin d) → ℝ => (Y i j) ^ m)
      (Measure.pi fun _ : Fin n => Measure.pi fun _ : Fin d => gaussianReal 0 1) :=
    PiLaw.integrable_comp_coord _ h1 i
  exact h2

/-! ### The operator norm: the trace bound and the entry bound -/

/-- A diagonal entry is at most the `l2` operator norm. -/
private theorem entry_le_l2_opNorm {D : ℕ} (M : Matrix (Fin D) (Fin D) ℝ) (i : Fin D) :
    M i i ≤ ‖M‖ := by
  classical
  set e : Fin D → ℝ := fun r => if r = i then (1 : ℝ) else 0 with he
  set x : EuclideanSpace ℝ (Fin D) := WithLp.toLp 2 e with hx
  have hxn : ‖x‖ = 1 := by
    have h1 : ‖x‖ ^ 2 = ∑ r, e r ^ 2 := EuclideanSpace.real_norm_sq_eq x
    have h2 : ∑ r, e r ^ 2 = 1 := by
      rw [Finset.sum_eq_single i] <;> simp +contextual [he]
    have h3 : ‖x‖ ^ 2 = 1 := by rw [h1, h2]
    nlinarith [norm_nonneg x]
  have hmv : M *ᵥ e = fun s => M s i := by
    funext s
    simp [Matrix.mulVec, dotProduct, he]
  have h := Matrix.l2_opNorm_mulVec M x
  rw [hxn, mul_one] at h
  have hval : ‖(EuclideanSpace.equiv (Fin D) ℝ).symm (M *ᵥ x.ofLp)‖ ^ 2 = ∑ s, (M s i) ^ 2 := by
    have hh : (M *ᵥ (x.ofLp : Fin D → ℝ)) = fun s => M s i := hmv
    rw [EuclideanSpace.real_norm_sq_eq]
    simp [hh]
  have hsingle : (M i i) ^ 2 ≤ ∑ s, (M s i) ^ 2 :=
    Finset.single_le_sum (f := fun s => (M s i) ^ 2) (fun s _ => sq_nonneg _) (Finset.mem_univ i)
  nlinarith [h, hval, hsingle, norm_nonneg ((EuclideanSpace.equiv (Fin D) ℝ).symm (M *ᵥ x.ofLp)),
    norm_nonneg M, abs_nonneg (M i i), sq_abs (M i i), le_abs_self (M i i)]

/-- `trace ((YᵀY)^k) ≤ d ‖Y‖^{2k}`: `d` diagonal entries, each at most `‖(YᵀY)^k‖`, and the
`l2` operator norm is submultiplicative with `‖YᵀY‖ = ‖Y‖²`. -/
private theorem trace_gram_pow_le {n d k : ℕ} (hk : 1 ≤ k) (Y : Matrix (Fin n) (Fin d) ℝ) :
    Matrix.trace ((Yᵀ * Y) ^ k) ≤ (d : ℝ) * ‖Y‖ ^ (2 * k) := by
  have h0 : Matrix.trace ((Yᵀ * Y) ^ k) = ∑ i, ((Yᵀ * Y) ^ k) i i := rfl
  have h1 : ∑ i, ((Yᵀ * Y) ^ k) i i ≤ ∑ _i : Fin d, ‖(Yᵀ * Y) ^ k‖ :=
    Finset.sum_le_sum fun i _ => entry_le_l2_opNorm _ i
  have h2 : ‖(Yᵀ * Y) ^ k‖ ≤ ‖Yᵀ * Y‖ ^ k := norm_pow_le' _ hk
  have h3 : ‖Yᵀ * Y‖ = ‖Y‖ * ‖Y‖ := by
    simpa using Matrix.l2_opNorm_conjTranspose_mul_self Y
  have h4 : (‖Y‖ * ‖Y‖) ^ k = ‖Y‖ ^ (2 * k) := by
    rw [← sq, ← pow_mul]
  have h5 : ∑ _i : Fin d, ‖(Yᵀ * Y) ^ k‖ = (d : ℝ) * ‖(Yᵀ * Y) ^ k‖ := by
    simp [Finset.sum_const, nsmul_eq_mul]
  rw [h0]
  calc ∑ i, ((Yᵀ * Y) ^ k) i i ≤ ∑ _i : Fin d, ‖(Yᵀ * Y) ^ k‖ := h1
    _ = (d : ℝ) * ‖(Yᵀ * Y) ^ k‖ := h5
    _ ≤ (d : ℝ) * ‖Y‖ ^ (2 * k) := by
        refine mul_le_mul_of_nonneg_left ?_ (Nat.cast_nonneg d)
        rw [← h4, ← h3]; exact h2

/-- `‖Y‖^{2(M+1)} ≤ (nd)^M ∑ entries^{2(M+1)}`: the operator norm is at most the Frobenius
norm (`R3.l2_opNorm_le_frobenius`), and Jensen bounds a power of a sum of squares. -/
private theorem opNorm_pow_le_sum {n d : ℕ} (Y : Matrix (Fin n) (Fin d) ℝ) (M : ℕ) :
    ‖Y‖ ^ (2 * (M + 1))
      ≤ ((n : ℝ) * d) ^ M * ∑ q : Fin n × Fin d, (Y q.1 q.2) ^ (2 * (M + 1)) := by
  classical
  set F : ℝ := ∑ q : Fin n × Fin d, (Y q.1 q.2) ^ 2 with hF
  have hF0 : 0 ≤ F := Finset.sum_nonneg fun q _ => sq_nonneg _
  have hFsum : ∑ i, ∑ j, Y i j ^ 2 = F := by
    rw [hF, Fintype.sum_prod_type]
  have hfrob : ‖Y‖ ≤ Real.sqrt F := by
    have h := R3.l2_opNorm_le_frobenius Y
    rwa [hFsum] at h
  have hsq : ‖Y‖ ^ 2 ≤ F := by
    have h := pow_le_pow_left₀ (norm_nonneg Y) hfrob 2
    rwa [Real.sq_sqrt hF0] at h
  have hstep1 : ‖Y‖ ^ (2 * (M + 1)) ≤ F ^ (M + 1) := by
    rw [pow_mul]
    exact pow_le_pow_left₀ (by positivity) hsq (M + 1)
  have hjen : F ^ (M + 1)
      ≤ ((Finset.univ : Finset (Fin n × Fin d)).card : ℝ) ^ M
        * ∑ q : Fin n × Fin d, ((Y q.1 q.2) ^ 2) ^ (M + 1) :=
    pow_sum_le_card_mul_sum_pow (fun q _ => sq_nonneg _) M
  have hcard : ((Finset.univ : Finset (Fin n × Fin d)).card : ℝ) = (n : ℝ) * d := by
    simp
  have hpow : ∀ q : Fin n × Fin d, ((Y q.1 q.2) ^ 2) ^ (M + 1) = (Y q.1 q.2) ^ (2 * (M + 1)) :=
    fun q => by rw [← pow_mul]
  rw [hcard] at hjen
  simp only [hpow] at hjen
  linarith

/-- Every even power of the operator norm is integrable under the Gaussian matrix law. -/
private theorem integrable_opNorm_pow (n d : ℕ) {k : ℕ} (hk : 1 ≤ k) :
    Integrable (fun Y : Matrix (Fin n) (Fin d) ℝ => ‖Y‖ ^ (2 * k)) (gaussianMatrix n d) := by
  obtain ⟨M, rfl⟩ : ∃ M, k = M + 1 := ⟨k - 1, by omega⟩
  have hg : Integrable (fun Y : Matrix (Fin n) (Fin d) ℝ =>
      ((n : ℝ) * d) ^ M * ∑ q : Fin n × Fin d, (Y q.1 q.2) ^ (2 * (M + 1)))
      (gaussianMatrix n d) :=
    (integrable_finsetSum _ fun q _ => integrable_entry_pow n d _ q.1 q.2).const_mul _
  refine Integrable.mono' hg ((R3.measurable_opNorm n d).pow_const _).aestronglyMeasurable ?_
  filter_upwards with Y
  rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
  exact opNorm_pow_le_sum Y M

/-- `∫ ‖Y‖^{2k} ≤ (nd)^k 4 (2k)!` under the Gaussian matrix law. -/
private theorem integral_opNorm_pow_le (n d : ℕ) {k : ℕ} (hk : 1 ≤ k) :
    ∫ Y, ‖Y‖ ^ (2 * k) ∂(gaussianMatrix n d)
      ≤ ((n : ℝ) * d) ^ k * (4 * (Nat.factorial (2 * k) : ℝ)) := by
  obtain ⟨M, rfl⟩ : ∃ M, k = M + 1 := ⟨k - 1, by omega⟩
  have hgi : Integrable (fun Y : Matrix (Fin n) (Fin d) ℝ =>
      ((n : ℝ) * d) ^ M * ∑ q : Fin n × Fin d, (Y q.1 q.2) ^ (2 * (M + 1)))
      (gaussianMatrix n d) :=
    (integrable_finsetSum _ fun q _ => integrable_entry_pow n d _ q.1 q.2).const_mul _
  have h1 : ∫ Y, ‖Y‖ ^ (2 * (M + 1)) ∂(gaussianMatrix n d)
      ≤ ∫ Y, ((n : ℝ) * d) ^ M * ∑ q : Fin n × Fin d, (Y q.1 q.2) ^ (2 * (M + 1))
          ∂(gaussianMatrix n d) :=
    integral_mono_of_nonneg (Eventually.of_forall fun Y => by positivity) hgi
      (Eventually.of_forall fun Y => opNorm_pow_le_sum Y M)
  have h2 : ∫ Y, ((n : ℝ) * d) ^ M * ∑ q : Fin n × Fin d, (Y q.1 q.2) ^ (2 * (M + 1))
        ∂(gaussianMatrix n d)
      = ((n : ℝ) * d) ^ M
        * ∑ q : Fin n × Fin d, ∫ Y, (Y q.1 q.2) ^ (2 * (M + 1)) ∂(gaussianMatrix n d) := by
    rw [integral_const_mul,
      integral_finsetSum _ fun q _ => integrable_entry_pow n d _ q.1 q.2]
  have h3 : ∑ q : Fin n × Fin d, ∫ Y, (Y q.1 q.2) ^ (2 * (M + 1)) ∂(gaussianMatrix n d)
      = ((n : ℝ) * d) * ∫ x : ℝ, x ^ (2 * (M + 1)) ∂(gaussianReal 0 1) := by
    rw [Finset.sum_congr rfl fun q _ => integral_entry_pow n d _ q.1 q.2]
    simp [Finset.sum_const, nsmul_eq_mul]
  have h4 : ∫ x : ℝ, x ^ (2 * (M + 1)) ∂(gaussianReal 0 1)
      ≤ 4 * (Nat.factorial (2 * (M + 1)) : ℝ) := integral_pow_gaussianReal_le _
  have hnn : (0 : ℝ) ≤ ((n : ℝ) * d) ^ M * ((n : ℝ) * d) := by positivity
  calc ∫ Y, ‖Y‖ ^ (2 * (M + 1)) ∂(gaussianMatrix n d)
      ≤ ((n : ℝ) * d) ^ M * (((n : ℝ) * d) * ∫ x : ℝ, x ^ (2 * (M + 1)) ∂(gaussianReal 0 1)) := by
        rw [h2, h3] at h1; exact h1
    _ = (((n : ℝ) * d) ^ M * ((n : ℝ) * d)) * ∫ x : ℝ, x ^ (2 * (M + 1)) ∂(gaussianReal 0 1) := by
        ring
    _ ≤ (((n : ℝ) * d) ^ M * ((n : ℝ) * d)) * (4 * (Nat.factorial (2 * (M + 1)) : ℝ)) :=
        mul_le_mul_of_nonneg_left h4 hnn
    _ = ((n : ℝ) * d) ^ (M + 1) * (4 * (Nat.factorial (2 * (M + 1)) : ℝ)) := by
        rw [pow_succ]

/-! ### The two scalar steps -/

/-- The split at `a`, with an arithmetic-geometric mean bound above `a`: below `a` the first
term wins, above `a` the value `ind ≥ 1` of the indicator pays for
`x^{2k} ≤ (θ/2) x^{4k} + 1/(2θ)`. -/
private theorem pow_le_split {x a θ ind : ℝ} (hx : 0 ≤ x) (ha : 0 ≤ a) (hθ : 0 < θ) (k : ℕ)
    (hind0 : 0 ≤ ind) (hind : a ≤ x → 1 ≤ ind) :
    x ^ (2 * k) ≤ a ^ (2 * k) + (θ / 2) * x ^ (2 * (2 * k)) + (1 / (2 * θ)) * ind := by
  have hxk : 0 ≤ x ^ (2 * k) := by positivity
  have hak : 0 ≤ a ^ (2 * k) := by positivity
  have hsq : x ^ (2 * (2 * k)) = (x ^ (2 * k)) ^ 2 := by
    rw [show 2 * (2 * k) = (2 * k) * 2 from by ring, pow_mul]
  rcases le_or_gt a x with h | h
  · have h1 : 1 ≤ ind := hind h
    have hid : (θ / 2) * (x ^ (2 * k)) ^ 2 + 1 / (2 * θ) - x ^ (2 * k)
        = (θ * x ^ (2 * k) - 1) ^ 2 / (2 * θ) := by
      field_simp
      ring
    have hnn : 0 ≤ (θ * x ^ (2 * k) - 1) ^ 2 / (2 * θ) := by positivity
    have h3 : (0 : ℝ) < 1 / (2 * θ) := by positivity
    have h2 : 1 / (2 * θ) ≤ (1 / (2 * θ)) * ind := by nlinarith [h3, h1]
    rw [hsq]
    linarith
  · have h1 : x ^ (2 * k) ≤ a ^ (2 * k) := pow_le_pow_left₀ hx h.le _
    have h2 : 0 ≤ (θ / 2) * x ^ (2 * (2 * k)) := by positivity
    have h3 : 0 ≤ (1 / (2 * θ)) * ind := by positivity
    linarith

/-- `log 4 + 10 k log D ≤ κ² D / 8` once `√D` passes the threshold `8 A / κ²` with
`A = log 4 + 160 C + 40`: at `k ≤ C log D + 1` the left side is `O(log² D)`, and
`log D ≤ 4 D^{1/4}` turns that into `O(√D)`. -/
private theorem log_bound {C κ D kk : ℝ} (hC : 0 < C) (hκ : 0 < κ) (hD : 1 ≤ D)
    (hk : kk ≤ C * Real.log D + 1)
    (hs : 8 * (Real.log 4 + 160 * C + 40) / κ ^ 2 ≤ Real.sqrt D) :
    Real.log 4 + 10 * kk * Real.log D ≤ κ ^ 2 * D / 8 := by
  have hD0 : (0 : ℝ) < D := lt_of_lt_of_le one_pos hD
  have hsqD0 : (0 : ℝ) ≤ Real.sqrt D := Real.sqrt_nonneg D
  set s : ℝ := Real.sqrt (Real.sqrt D) with hsdef
  have hs0 : 0 ≤ s := Real.sqrt_nonneg _
  have hs2 : s ^ 2 = Real.sqrt D := Real.sq_sqrt hsqD0
  have hs4 : s ^ 4 = D := by
    have h1 : (s ^ 2) ^ 2 = D := by rw [hs2, Real.sq_sqrt hD0.le]
    nlinarith [h1]
  have hs1 : 1 ≤ s := by
    have h1 : (1 : ℝ) ≤ Real.sqrt D := by
      have h := Real.sqrt_le_sqrt hD
      simpa using h
    have h2 := Real.sqrt_le_sqrt h1
    simpa [hsdef] using h2
  have hL0 : 0 ≤ Real.log D := Real.log_nonneg hD
  have hLs : Real.log D ≤ 4 * s := by
    have h1 : Real.log D = 4 * Real.log s := by
      rw [← hs4, Real.log_pow]; push_cast; ring
    have h2 : Real.log s ≤ s - 1 := Real.log_le_sub_one_of_pos (by linarith)
    linarith
  have hlog4 : (0 : ℝ) ≤ Real.log 4 := Real.log_nonneg (by norm_num)
  have hL2 : Real.log D ^ 2 ≤ 16 * s ^ 2 := by nlinarith [hLs, hL0]
  have hss : s ≤ s ^ 2 := by nlinarith [hs1]
  have h1s : (1 : ℝ) ≤ s ^ 2 := by nlinarith [hs1]
  have e1 : 10 * kk * Real.log D ≤ 10 * C * Real.log D ^ 2 + 10 * Real.log D := by
    nlinarith [mul_le_mul_of_nonneg_left hk hL0]
  have e2 : 10 * C * Real.log D ^ 2 ≤ 160 * C * s ^ 2 := by nlinarith [hL2, hC.le]
  have e3 : 10 * Real.log D ≤ 40 * s ^ 2 := by nlinarith [hLs, hss, hs0]
  have e4 : Real.log 4 ≤ Real.log 4 * s ^ 2 := by nlinarith [h1s, hlog4]
  have hstep : Real.log 4 + 10 * kk * Real.log D
      ≤ (Real.log 4 + 160 * C + 40) * s ^ 2 := by nlinarith [e1, e2, e3, e4]
  have hthr : 8 * (Real.log 4 + 160 * C + 40) ≤ κ ^ 2 * s ^ 2 := by
    rw [div_le_iff₀ (by positivity : (0 : ℝ) < κ ^ 2), ← hs2] at hs
    linarith
  have hfin : (Real.log 4 + 160 * C + 40) * s ^ 2 ≤ κ ^ 2 * D / 8 := by
    rw [← hs4]
    nlinarith [mul_le_mul_of_nonneg_right hthr (sq_nonneg s)]
  linarith

/-! ### The bound -/

set_option maxHeartbeats 1000000 in
-- One declaration assembles nine estimates (the split, the fourth moment, the tail and the
-- trace step), so the default 200000 heartbeats is not enough.
/-- G. For `N` large, at `k = momOrder C (d N)`,
`∫ trace ((Y_GᵀY_G)^k) ≤ 2 d (d (bulkEdge c + ε))^k`. -/
theorem gaussian_trace_le {c : ℝ} (hc : 0 < c) {n d : ℕ → ℕ}
    (hn : ∀ N, 0 < n N) (hd : ∀ N, 0 < d N) (hdtop : Tendsto d atTop atTop)
    (hcN : Tendsto (fun N => (n N : ℝ) / d N) atTop (𝓝 c)) {C : ℝ} (hC : 0 < C)
    {ε : ℝ} (hε : 0 < ε) :
    ∀ᶠ N in atTop,
      ∫ Y, Matrix.trace ((Yᵀ * Y) ^ (momOrder C (d N))) ∂(gaussianMatrix (n N) (d N))
        ≤ 2 * (d N : ℝ) * ((d N : ℝ) * (bulkEdge c + ε)) ^ (momOrder C (d N)) := by
  classical
  -- the concentration scale `κ`, exactly as in `R3.tendsto_measure_lamMax_le`
  set sc : ℝ := Real.sqrt c with hscdef
  have hsc0 : (0 : ℝ) ≤ sc := by rw [hscdef]; exact Real.sqrt_nonneg c
  have hbe : bulkEdge c = (1 + sc) ^ 2 := by rw [hscdef]; rfl
  clear_value sc
  have hK0 : (0 : ℝ) < 4 * (1 + sc) + 4 := by nlinarith
  set κ : ℝ := min 1 (ε / (4 * (1 + sc) + 4)) with hkapdef
  have hκ0 : (0 : ℝ) < κ := by rw [hkapdef]; exact lt_min one_pos (div_pos hε hK0)
  have hκ1 : κ ≤ 1 := by rw [hkapdef]; exact min_le_left _ _
  have hκK : κ * (4 * (1 + sc) + 4) ≤ ε := by
    have hle : κ ≤ ε / (4 * (1 + sc) + 4) := by rw [hkapdef]; exact min_le_right _ _
    have h := mul_le_mul_of_nonneg_right hle hK0.le
    rwa [div_mul_cancel₀ _ hK0.ne'] at h
  clear_value κ
  have hlog4 : (0 : ℝ) ≤ Real.log 4 := Real.log_nonneg (by norm_num)
  have hAnn : (0 : ℝ) ≤ 8 * (Real.log 4 + 160 * C + 40) / κ ^ 2 :=
    div_nonneg (by nlinarith [hC.le]) (by positivity)
  -- the eventual hypotheses
  have hdR : Tendsto (fun N => (d N : ℝ)) atTop atTop := tendsto_natCast_atTop_atTop.comp hdtop
  have E1 := eventually_one_le_momOrder (d := d) hdtop hC
  have E2 := eventually_momOrder_le (d := d) hdtop hC
  have E3 : ∀ᶠ N in atTop, Real.sqrt ((n N : ℝ) / (d N : ℝ)) < sc + κ :=
    Filter.Tendsto.eventually_lt_const (by linarith) hcN.sqrt
  have E4 : ∀ᶠ N in atTop, (n N : ℝ) / (d N : ℝ) < c + 1 :=
    Filter.Tendsto.eventually_lt_const (by linarith) hcN
  have E5 : ∀ᶠ N in atTop, (16 * (c + 1) : ℝ) ≤ (d N : ℝ) := hdR.eventually_ge_atTop _
  have E6 : ∀ᶠ N in atTop, (1 : ℝ) ≤ (d N : ℝ) := hdR.eventually_ge_atTop _
  have E7 : ∀ᶠ N in atTop,
      ((8 * (Real.log 4 + 160 * C + 40) / κ ^ 2) ^ 2 : ℝ) ≤ (d N : ℝ) := hdR.eventually_ge_atTop _
  filter_upwards [E1, E2, E3, E4, E5, E6, E7] with N hk1 hk8 hrN hnN h16 hd1 hA7
  have hdd0 : (0 : ℝ) < (d N : ℝ) := by exact_mod_cast hd N
  have hnn0 : (0 : ℝ) < (n N : ℝ) := by exact_mod_cast hn N
  -- the moment order, made opaque once its two arithmetic facts are out
  set k : ℕ := momOrder C (d N) with hkdef
  have hkd : (k : ℝ) ≤ (d N : ℝ) := by
    have h : k ≤ d N := by omega
    exact_mod_cast h
  have hklog : (k : ℝ) ≤ C * Real.log (d N : ℝ) + 1 := by
    have hpos : (0 : ℝ) ≤ C * Real.log (d N : ℝ) := mul_nonneg hC.le (Real.log_nonneg hd1)
    have h := Nat.ceil_lt_add_one (a := C * Real.log (d N : ℝ)) hpos
    simp only [hkdef, momOrder]
    exact h.le
  clear_value k
  have hk0R : (0 : ℝ) ≤ (k : ℝ) := Nat.cast_nonneg k
  -- the two facts of item R3 about the Gaussian operator norm
  have hsd0 : (0 : ℝ) < Real.sqrt (d N : ℝ) := Real.sqrt_pos.mpr hdd0
  have hsdsq : Real.sqrt (d N : ℝ) ^ 2 = (d N : ℝ) := Real.sq_sqrt hdd0.le
  have hr0 : (0 : ℝ) ≤ Real.sqrt ((n N : ℝ) / (d N : ℝ)) := Real.sqrt_nonneg _
  have hsp : Real.sqrt ((n N : ℝ) / (d N : ℝ)) * Real.sqrt (d N : ℝ) = Real.sqrt (n N : ℝ) := by
    rw [Real.sqrt_div (le_of_lt hnn0), div_mul_cancel₀ _ hsd0.ne']
  have hGordon := R3.integral_opNorm_le (p := n N) (d := d N) (hn N) (hd N)
  have hStail := R3.measure_opNorm_ge_le (p := n N) (d := d N) (hn N) (hd N)
    (show (0 : ℝ) < κ * Real.sqrt (d N : ℝ) from mul_pos hκ0 hsd0)
  have hμ0 : (0 : ℝ) ≤ ∫ Z', ‖Z'‖ ∂(gaussianMatrix (n N) (d N)) :=
    integral_nonneg fun Z => norm_nonneg _
  -- the truncation level `a`, made opaque once its two facts are out
  set a : ℝ := (∫ Z', ‖Z'‖ ∂(gaussianMatrix (n N) (d N))) + κ * Real.sqrt (d N : ℝ) with hadef
  have ha0 : 0 ≤ a := by
    rw [hadef]; exact add_nonneg hμ0 (mul_nonneg hκ0.le hsd0.le)
  have hale : a ≤ (Real.sqrt ((n N : ℝ) / (d N : ℝ)) + 1 + κ) * Real.sqrt (d N : ℝ) := by
    rw [hadef]; nlinarith [hGordon, hsp]
  clear_value a
  clear hadef hμ0 hGordon
  -- the bad event, made opaque once its two facts are out
  set S : Set (Matrix (Fin (n N)) (Fin (d N)) ℝ) := {Z | a ≤ ‖Z‖} with hSdef
  have hSmeas : MeasurableSet S := by
    rw [hSdef]; exact measurableSet_le measurable_const (R3.measurable_opNorm _ _)
  have hSmem : ∀ Z : Matrix (Fin (n N)) (Fin (d N)) ℝ, a ≤ ‖Z‖ → Z ∈ S := by
    intro Z hZ; rw [hSdef]; exact hZ
  have hStail' : ((gaussianMatrix (n N) (d N)) S).toReal
      ≤ Real.exp (-(κ ^ 2 * (d N : ℝ)) / 2) := by
    have hsq : (κ * Real.sqrt (d N : ℝ)) ^ 2 = κ ^ 2 * (d N : ℝ) := by rw [mul_pow, hsdsq]
    rwa [hsq] at hStail
  clear_value S
  clear hSdef hStail
  -- the tail scale `θ`, made opaque once its three facts are out
  set θ : ℝ := Real.exp (-(κ ^ 2 * (d N : ℝ)) / 8) with hθdef
  have hθ0 : 0 < θ := by rw [hθdef]; exact Real.exp_pos _
  have hE : Real.exp ((κ ^ 2 * (d N : ℝ)) / 8) * θ = 1 := by
    have hz : (κ ^ 2 * (d N : ℝ)) / 8 + -(κ ^ 2 * (d N : ℝ)) / 8 = 0 := by ring
    rw [hθdef, ← Real.exp_add, hz, Real.exp_zero]
  have hθtail : Real.exp (-(κ ^ 2 * (d N : ℝ)) / 2) ≤ θ := by
    rw [hθdef]
    exact Real.exp_le_exp.mpr (by nlinarith [sq_nonneg κ, hdd0])
  clear_value θ
  clear hθdef
  -- the target constant
  have hbe1 : (1 : ℝ) ≤ bulkEdge c + ε := by rw [hbe]; nlinarith [hsc0, hε]
  have hP1 : (1 : ℝ) ≤ ((d N : ℝ) * (bulkEdge c + ε)) ^ k :=
    one_le_pow₀ (by nlinarith [hd1, hbe1])
  -- step 1: the truncation level is below the edge
  have hterm1 : a ^ (2 * k) ≤ ((d N : ℝ) * (bulkEdge c + ε)) ^ k := by
    have hedge := R3.edge_arith hsc0 hκ0 hκ1 hκK hr0 hrN.le
    have hasq : a ^ 2 ≤ (d N : ℝ) * (bulkEdge c + ε) := by
      have h1 : a ^ 2
          ≤ ((Real.sqrt ((n N : ℝ) / (d N : ℝ)) + 1 + κ) * Real.sqrt (d N : ℝ)) ^ 2 :=
        pow_le_pow_left₀ ha0 hale 2
      have h2 : ((Real.sqrt ((n N : ℝ) / (d N : ℝ)) + 1 + κ) * Real.sqrt (d N : ℝ)) ^ 2
          = (Real.sqrt ((n N : ℝ) / (d N : ℝ)) + 1 + κ) ^ 2 * (d N : ℝ) := by
        rw [mul_pow, hsdsq]
      have h3 : (Real.sqrt ((n N : ℝ) / (d N : ℝ)) + 1 + κ) ^ 2 * (d N : ℝ)
          ≤ ((1 + sc) ^ 2 + ε) * (d N : ℝ) := mul_le_mul_of_nonneg_right hedge hdd0.le
      rw [hbe]
      linarith
    calc a ^ (2 * k) = (a ^ 2) ^ k := by rw [pow_mul]
      _ ≤ ((d N : ℝ) * (bulkEdge c + ε)) ^ k := pow_le_pow_left₀ (sq_nonneg a) hasq k
  -- step 2: the fourth moment against the tail scale
  have hterm2 : (θ / 2) * (∫ Y, ‖Y‖ ^ (2 * (2 * k)) ∂(gaussianMatrix (n N) (d N))) ≤ 1 / 2 := by
    have hI := integral_opNorm_pow_le (n N) (d N) (k := 2 * k) (by omega)
    have hfacle : (Nat.factorial (2 * (2 * k)) : ℝ) ≤ (16 * (k : ℝ) ^ 2) ^ (2 * k) := by
      have h1 : Nat.factorial (2 * (2 * k)) ≤ (2 * (2 * k)) ^ (2 * (2 * k)) :=
        Nat.factorial_le_pow _
      have h1' : (Nat.factorial (2 * (2 * k)) : ℝ)
          ≤ (((2 * (2 * k) : ℕ) : ℝ)) ^ (2 * (2 * k)) := by exact_mod_cast h1
      have h2 : (((2 * (2 * k) : ℕ) : ℝ)) ^ (2 * (2 * k)) = (16 * (k : ℝ) ^ 2) ^ (2 * k) := by
        have hcast : ((2 * (2 * k) : ℕ) : ℝ) = 4 * (k : ℝ) := by push_cast; ring
        rw [hcast, pow_mul]
        congr 1
        ring
      rwa [h2] at h1'
    have hnd : (n N : ℝ) ≤ (c + 1) * (d N : ℝ) := by
      rw [div_lt_iff₀ hdd0] at hnN
      linarith
    have hk2 : (k : ℝ) ^ 2 ≤ (d N : ℝ) ^ 2 := pow_le_pow_left₀ hk0R hkd 2
    have hbig : ((n N : ℝ) * (d N : ℝ)) * (16 * (k : ℝ) ^ 2) ≤ ((d N : ℝ)) ^ 5 := by
      have s1 : ((n N : ℝ) * (d N : ℝ)) * (16 * (k : ℝ) ^ 2)
          ≤ (((c + 1) * (d N : ℝ)) * (d N : ℝ)) * (16 * (d N : ℝ) ^ 2) := by
        refine mul_le_mul (mul_le_mul_of_nonneg_right hnd hdd0.le) (by linarith) (by positivity)
          (mul_nonneg (mul_nonneg (by linarith) hdd0.le) hdd0.le)
      have s2 : (((c + 1) * (d N : ℝ)) * (d N : ℝ)) * (16 * (d N : ℝ) ^ 2)
          = (16 * (c + 1)) * (d N : ℝ) ^ 4 := by ring
      have s3 : (16 * (c + 1)) * (d N : ℝ) ^ 4 ≤ (d N : ℝ) * (d N : ℝ) ^ 4 :=
        mul_le_mul_of_nonneg_right h16 (by positivity)
      have s4 : (d N : ℝ) * (d N : ℝ) ^ 4 = ((d N : ℝ)) ^ 5 := by ring
      linarith
    have hB : ((n N : ℝ) * (d N : ℝ)) ^ (2 * k) * (4 * (Nat.factorial (2 * (2 * k)) : ℝ))
        ≤ 4 * (((d N : ℝ)) ^ 5) ^ (2 * k) := by
      have hstep1 : ((n N : ℝ) * (d N : ℝ)) ^ (2 * k) * (4 * (Nat.factorial (2 * (2 * k)) : ℝ))
          ≤ ((n N : ℝ) * (d N : ℝ)) ^ (2 * k) * (4 * (16 * (k : ℝ) ^ 2) ^ (2 * k)) :=
        mul_le_mul_of_nonneg_left (by linarith) (by positivity)
      have hstep2 : ((n N : ℝ) * (d N : ℝ)) ^ (2 * k) * (4 * (16 * (k : ℝ) ^ 2) ^ (2 * k))
          = 4 * (((n N : ℝ) * (d N : ℝ)) * (16 * (k : ℝ) ^ 2)) ^ (2 * k) := by
        rw [mul_pow]; ring
      have hstep3 : (((n N : ℝ) * (d N : ℝ)) * (16 * (k : ℝ) ^ 2)) ^ (2 * k)
          ≤ (((d N : ℝ)) ^ 5) ^ (2 * k) := pow_le_pow_left₀ (by positivity) hbig _
      linarith
    have hA7' : 8 * (Real.log 4 + 160 * C + 40) / κ ^ 2 ≤ Real.sqrt (d N : ℝ) := by
      have h := Real.sqrt_le_sqrt hA7
      rwa [Real.sqrt_sq hAnn] at h
    have hlog := log_bound hC hκ0 hd1 hklog hA7'
    have hexpeq : 4 * (((d N : ℝ)) ^ 5) ^ (2 * k)
        = Real.exp (Real.log 4 + 10 * (k : ℝ) * Real.log (d N : ℝ)) := by
      rw [Real.exp_add, Real.exp_log (by norm_num : (0 : ℝ) < 4)]
      congr 1
      have hrw : 10 * (k : ℝ) * Real.log (d N : ℝ)
          = ((2 * k : ℕ) : ℝ) * Real.log ((d N : ℝ) ^ 5) := by
        rw [Real.log_pow]; push_cast; ring
      rw [hrw, Real.exp_nat_mul, Real.exp_log (by positivity)]
    have hBE : ((n N : ℝ) * (d N : ℝ)) ^ (2 * k) * (4 * (Nat.factorial (2 * (2 * k)) : ℝ))
        ≤ Real.exp ((κ ^ 2 * (d N : ℝ)) / 8) := by
      refine hB.trans ?_
      rw [hexpeq]
      exact Real.exp_le_exp.mpr (by linarith)
    have hfinal : (θ / 2) * (∫ Y, ‖Y‖ ^ (2 * (2 * k)) ∂(gaussianMatrix (n N) (d N)))
        ≤ (θ / 2) * Real.exp ((κ ^ 2 * (d N : ℝ)) / 8) :=
      mul_le_mul_of_nonneg_left (hI.trans hBE) (by positivity)
    have hval : (θ / 2) * Real.exp ((κ ^ 2 * (d N : ℝ)) / 8) = 1 / 2 := by
      rw [show (θ / 2) * Real.exp ((κ ^ 2 * (d N : ℝ)) / 8)
          = (Real.exp ((κ ^ 2 * (d N : ℝ)) / 8) * θ) / 2 from by ring, hE]
    linarith
  -- step 3: the indicator term
  have hterm3 : (1 / (2 * θ)) * ((gaussianMatrix (n N) (d N)) S).toReal ≤ 1 / 2 := by
    have hpos : (0 : ℝ) ≤ 1 / (2 * θ) := by positivity
    have h1 : (1 / (2 * θ)) * ((gaussianMatrix (n N) (d N)) S).toReal ≤ (1 / (2 * θ)) * θ :=
      mul_le_mul_of_nonneg_left (hStail'.trans hθtail) hpos
    have h2 : (1 / (2 * θ)) * θ = 1 / 2 := by field_simp
    linarith
  -- the split, integrated
  have hint2 : Integrable (fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ => ‖Y‖ ^ (2 * (2 * k)))
      (gaussianMatrix (n N) (d N)) := integrable_opNorm_pow _ _ (by omega)
  have hintInd : Integrable
      (Set.indicator S (1 : Matrix (Fin (n N)) (Fin (d N)) ℝ → ℝ))
      (gaussianMatrix (n N) (d N)) := (integrable_const (1 : ℝ)).indicator hSmeas
  have hi1 : Integrable (fun _Y : Matrix (Fin (n N)) (Fin (d N)) ℝ => a ^ (2 * k))
      (gaussianMatrix (n N) (d N)) := integrable_const _
  have hi2 : Integrable
      (fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ => (θ / 2) * ‖Y‖ ^ (2 * (2 * k)))
      (gaussianMatrix (n N) (d N)) := hint2.const_mul _
  have hi3 : Integrable (fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
      (1 / (2 * θ)) * Set.indicator S (1 : Matrix (Fin (n N)) (Fin (d N)) ℝ → ℝ) Y)
      (gaussianMatrix (n N) (d N)) := hintInd.const_mul _
  have hi12 : Integrable (fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
      a ^ (2 * k) + (θ / 2) * ‖Y‖ ^ (2 * (2 * k))) (gaussianMatrix (n N) (d N)) := hi1.add hi2
  have hintRHS : Integrable (fun Y : Matrix (Fin (n N)) (Fin (d N)) ℝ =>
      a ^ (2 * k) + (θ / 2) * ‖Y‖ ^ (2 * (2 * k))
        + (1 / (2 * θ)) * Set.indicator S (1 : Matrix (Fin (n N)) (Fin (d N)) ℝ → ℝ) Y)
      (gaussianMatrix (n N) (d N)) := hi12.add hi3
  have hsplit : ∀ Y : Matrix (Fin (n N)) (Fin (d N)) ℝ,
      ‖Y‖ ^ (2 * k) ≤ a ^ (2 * k) + (θ / 2) * ‖Y‖ ^ (2 * (2 * k))
        + (1 / (2 * θ)) * Set.indicator S (1 : Matrix (Fin (n N)) (Fin (d N)) ℝ → ℝ) Y := by
    intro Y
    refine pow_le_split (norm_nonneg Y) ha0 hθ0 k ?_ ?_
    · exact Set.indicator_nonneg (fun x _ => zero_le_one) Y
    · intro hY
      simp [Set.indicator_of_mem (hSmem Y hY)]
  have hmono : ∫ Y, ‖Y‖ ^ (2 * k) ∂(gaussianMatrix (n N) (d N))
      ≤ ∫ Y, (a ^ (2 * k) + (θ / 2) * ‖Y‖ ^ (2 * (2 * k))
          + (1 / (2 * θ)) * Set.indicator S (1 : Matrix (Fin (n N)) (Fin (d N)) ℝ → ℝ) Y)
        ∂(gaussianMatrix (n N) (d N)) :=
    integral_mono_of_nonneg (Eventually.of_forall fun Y => by positivity) hintRHS
      (Eventually.of_forall hsplit)
  have hRHS : ∫ Y, (a ^ (2 * k) + (θ / 2) * ‖Y‖ ^ (2 * (2 * k))
        + (1 / (2 * θ)) * Set.indicator S (1 : Matrix (Fin (n N)) (Fin (d N)) ℝ → ℝ) Y)
      ∂(gaussianMatrix (n N) (d N))
      = a ^ (2 * k) + (θ / 2) * (∫ Y, ‖Y‖ ^ (2 * (2 * k)) ∂(gaussianMatrix (n N) (d N)))
        + (1 / (2 * θ)) * ((gaussianMatrix (n N) (d N)) S).toReal := by
    have e1 : ∫ Y, (a ^ (2 * k) + (θ / 2) * ‖Y‖ ^ (2 * (2 * k))
          + (1 / (2 * θ)) * Set.indicator S (1 : Matrix (Fin (n N)) (Fin (d N)) ℝ → ℝ) Y)
        ∂(gaussianMatrix (n N) (d N))
        = (∫ Y, (a ^ (2 * k) + (θ / 2) * ‖Y‖ ^ (2 * (2 * k))) ∂(gaussianMatrix (n N) (d N)))
          + ∫ Y, ((1 / (2 * θ))
              * Set.indicator S (1 : Matrix (Fin (n N)) (Fin (d N)) ℝ → ℝ) Y)
            ∂(gaussianMatrix (n N) (d N)) := integral_add hi12 hi3
    have e2 : ∫ Y, (a ^ (2 * k) + (θ / 2) * ‖Y‖ ^ (2 * (2 * k))) ∂(gaussianMatrix (n N) (d N))
        = (∫ _Y : Matrix (Fin (n N)) (Fin (d N)) ℝ, a ^ (2 * k) ∂(gaussianMatrix (n N) (d N)))
          + ∫ Y, ((θ / 2) * ‖Y‖ ^ (2 * (2 * k))) ∂(gaussianMatrix (n N) (d N)) :=
      integral_add hi1 hi2
    have e3 : ∫ Y, ((θ / 2) * ‖Y‖ ^ (2 * (2 * k))) ∂(gaussianMatrix (n N) (d N))
        = (θ / 2) * ∫ Y, ‖Y‖ ^ (2 * (2 * k)) ∂(gaussianMatrix (n N) (d N)) :=
      integral_const_mul _ _
    have e4 : ∫ Y, ((1 / (2 * θ))
          * Set.indicator S (1 : Matrix (Fin (n N)) (Fin (d N)) ℝ → ℝ) Y)
        ∂(gaussianMatrix (n N) (d N))
        = (1 / (2 * θ)) * ∫ Y, Set.indicator S (1 : Matrix (Fin (n N)) (Fin (d N)) ℝ → ℝ) Y
            ∂(gaussianMatrix (n N) (d N)) := integral_const_mul _ _
    have e5 : ∫ Y, Set.indicator S (1 : Matrix (Fin (n N)) (Fin (d N)) ℝ → ℝ) Y
          ∂(gaussianMatrix (n N) (d N))
        = ((gaussianMatrix (n N) (d N)) S).toReal := by
      rw [integral_indicator_one hSmeas, measureReal_def]
    have e6 : ∫ _Y : Matrix (Fin (n N)) (Fin (d N)) ℝ, a ^ (2 * k)
        ∂(gaussianMatrix (n N) (d N)) = a ^ (2 * k) := by simp
    rw [e5] at e4
    linarith [e1, e2, e3, e4, e6]
  have hopbound : ∫ Y, ‖Y‖ ^ (2 * k) ∂(gaussianMatrix (n N) (d N))
      ≤ 2 * ((d N : ℝ) * (bulkEdge c + ε)) ^ k := by
    rw [hRHS] at hmono
    linarith
  -- from the operator norm to the trace
  have htrace : ∫ Y, Matrix.trace ((Yᵀ * Y) ^ k) ∂(gaussianMatrix (n N) (d N))
      ≤ (d N : ℝ) * ∫ Y, ‖Y‖ ^ (2 * k) ∂(gaussianMatrix (n N) (d N)) := by
    have h := integral_mono_of_nonneg
      (Eventually.of_forall fun Y => trace_gram_pow_nonneg (k := k) Y)
      ((integrable_opNorm_pow (n N) (d N) hk1).const_mul (d N : ℝ))
      (Eventually.of_forall fun Y => trace_gram_pow_le hk1 Y)
    rwa [integral_const_mul] at h
  calc ∫ Y, Matrix.trace ((Yᵀ * Y) ^ k) ∂(gaussianMatrix (n N) (d N))
      ≤ (d N : ℝ) * ∫ Y, ‖Y‖ ^ (2 * k) ∂(gaussianMatrix (n N) (d N)) := htrace
    _ ≤ (d N : ℝ) * (2 * ((d N : ℝ) * (bulkEdge c + ε)) ^ k) :=
        mul_le_mul_of_nonneg_left hopbound hdd0.le
    _ = 2 * (d N : ℝ) * ((d N : ℝ) * (bulkEdge c + ε)) ^ k := by ring

/-- G0. The Gaussian trace moment is nonnegative: `(YᵀY)^k` is positive semidefinite. -/
theorem integral_trace_gaussian_nonneg (n d k : ℕ) :
    0 ≤ ∫ Y, Matrix.trace ((Yᵀ * Y) ^ k) ∂(gaussianMatrix n d) :=
  integral_nonneg fun Y => trace_gram_pow_nonneg Y

end Edge
end StackedSVD
