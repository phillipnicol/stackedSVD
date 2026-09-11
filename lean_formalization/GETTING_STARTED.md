# Verifying a result in Lean: a beginner's guide

This guide shows how to install Lean and confirm that one theorem from this repository
compiles (i.e., the proof checker accepts it with no errors).

No prior Lean experience is required.

Written by Phillip Nicol (2026-09-09, as `example.md` on the branch `pnicol_edits`); moved
here and updated on 2026-09-10 with the numbers of a fresh clone-and-build test on a Linux
server. The theorem text below is the source of
`lean/StackedSVD/Main.lean` at that date; if the line numbers drift, `grep -n
thm_svdstack_weighted lean/StackedSVD/Main.lean` finds it.

---

## 1. Install Lean

Lean is installed through **elan**, its version manager. Run this in your terminal:

```sh
curl https://elan.lean-lang.org/elan-init.sh -sSf | sh
```

Follow the prompts. When it finishes, restart your terminal (or run `source ~/.profile`)
so that `elan` and `lake` are on your PATH.

Optional but recommended: install the **VS Code extension** for Lean 4, which shows proof
goals and errors inline as you edit.

```sh
code --install-extension leanprover.lean4
```

---

## 2. Clone this repository and download the pre-built library cache

The proofs depend on **Mathlib**, Lean's standard library of mathematics. Building Mathlib
from source takes hours; the Mathlib project ships a pre-built binary cache, and the
command below downloads the cache of the exact Mathlib version this repository pins.

```sh
git clone https://github.com/phillipnicol/stackedSVD.git
cd stackedSVD/lean_formalization/lean
lake exe cache get      # downloads pre-built Mathlib binaries (about 1 GB, 7 GB once unpacked)
```

`lake` is Lean's build tool (analogous to `make` or `cargo`). The `cache get` step usually
takes a few minutes on a fast connection. The second dependency, StatsMLlib (Weyl,
Davis-Kahan and Gaussian concentration), has no binary cache; `lake build` compiles it from
its pinned commit on the way (41 modules, counted in the build time below).

---

## 3. Build the project

```sh
lake build
```

This compiles all 230 modules of the project (210 source files of ours and 20 vendored
ones). The first build takes a while. A fresh clone built in 13 min 21 s on 6 cores of a
Linux server on 2026-09-10 (8975 Lake jobs, before the three modules of the count proof
were added the same day; the tree of commit `416ab9d` has 8978 jobs); expect 20 to 30 minutes on a laptop with 4
cores. Subsequent builds are fast because Lake caches compiled modules, and a build with no
source change replays in under a minute. The compiled files take about 10 GB of disk under
`lean/.lake/` (7.1 GB of them the unpacked Mathlib cache).

When it finishes with no errors, every theorem in the project has been accepted by the
Lean kernel.

---

## 4. A worked example: the optimal weighted SVD-stack theorem

**What the theorem says in words.**
You have M data tables, each of size n_i × d, each containing a shared rank-one signal hidden
in Gaussian noise. Some tables have a strong enough signal to be detectable above the noise
floor (the BBP threshold); call those the "active" tables. The paper's main question is: what
is the best way to combine these M tables to recover the shared signal direction v?

The answer given by `thm:svdstack_weighted` (the second theorem of Section 4 of the paper,
`main_paper.tex` line 499) is: stack the tables with the optimal weights

```
w_i = θ_i √((θ_i² + 1) / (θ_i² + c_i))   if θ_i⁴ > c_i  (active table)
w_i = 0                                     otherwise
```

and compute the top right singular vector of the weighted stack. Then, as n_i and d grow
proportionally (n_i / d → c_i), the squared alignment of the estimator with the true v
converges in probability to

```
S / (S + 1),   where S = Σ_{active i} β_i² / (1 − β_i²)
```

Here β_i is the "effective signal-to-noise ratio" of table i after the random-matrix
correction. The formula aggregates the contribution of every active table. The companion
theorem `thm_svdstack_weighted_opt` (the same file, lines 252 to 265) says this limit is
optimal: no other weighting does better.

**The Lean statement.**
Here is the statement as it appears in `lean/StackedSVD/Main.lean` (lines 238 to 245),
with the binders spread over more lines and a comment added on each hypothesis:

```lean
theorem thm_svdstack_weighted {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    {μ : ∀ N, Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i)                          -- every aspect ratio is positive
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))  -- β_i is defined by the model
    (hthr : ∃ k, 0 < β k)                        -- at least one table is detectable
    (hreg : ∀ i, (m.tbl i).Regime (c i))          -- proportional-growth regime
    (hG : m.JointGaussianNoise) :                 -- independent Gaussian noise
    TendstoInProb μ
      (fun N ω => m.svdstackPerfW (optW β) N ω)  -- ⟨v̂, v⟩² with optimal weights
      (svdstackLimitOpt β)                        -- → S / (S + 1)
```

Reading the hypotheses: `hc` says every table has a positive aspect ratio c_i; `hβdef` ties
the free parameter β_i to the model's signal strength and aspect ratio; `hthr` says at least
one table exceeds the detection threshold; `hreg` says the sample size n_i(N) and dimension
d(N) both grow and their ratio converges to c_i; `hG` says all tables have i.i.d. N(0, 1)
noise, independent across tables. The conclusion `TendstoInProb` is convergence in
probability of the squared alignment to the optimal limit.

**Checking it yourself.**
After building the project (step 3), make sure you are inside the `lean/` directory.
The repository already holds a file `check_weighted.lean` there with the four lines below; to
make your own, note that the contents must go in a file, not typed directly into the
terminal:

```sh
cat > check_weighted.lean << 'EOF'
import StackedSVD.Main

-- Print the type of the theorem; Lean will error if it is not accepted.
#check @StackedSVD.Main.thm_svdstack_weighted
EOF
```

Then run:

```sh
lake env lean check_weighted.lean
```

Expected output (Lean prints the theorem's type signature and exits with no error; the
leading `@` is Lean's mark for a name with every implicit argument shown):

```
@StackedSVD.Main.thm_svdstack_weighted : ∀ {Ω : ℕ → Type u_1} [inst : (N : ℕ) → MeasurableSpace (Ω N)]
  {μ : (N : ℕ) → MeasureTheory.Measure (Ω N)} {M : ℕ} {n : Fin M → ℕ → ℕ} {d : ℕ → ℕ} [inst_1 : NeZero M]
  ...
```

If you see that output with no `error:` lines, the kernel has confirmed the proof. Loading
the imports takes 10 s to 30 s; the check itself is instant.

**Every theorem at once.** `lean/StackedSVD/Main.lean` is the statement file: 40
theorems, one per paper result, each proved by one call into the tree. `docs/THEOREMS.md`
quotes each one next to the paper's LaTeX statement and explains its hypotheses in words.
The eight gate scripts of `scripts/` (`scripts/README.md`) check, among other things, that
no proof is left open and that no axiom beyond Lean's three standard ones is used.

---

