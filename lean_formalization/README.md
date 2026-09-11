# Stacked SVD or SVD stacked: the Lean formalization

## Summary

This folder (`lean_formalization` of the `stackedSVD` repository) holds a machine-checked
proof, in Lean 4 with Mathlib, of the results of
*"Stacked SVD or SVD stacked? A random matrix theory perspective on data integration"*
(Tavor Z. Baharav, Phillip B. Nicol, Rafael A. Irizarry and Rong Ma; arXiv:2507.22170).
The Lean checker accepts the tree with 0 open proofs and with no axiom beyond the 3
standard ones of Lean. Every result of the paper holds under the hypotheses that its row of
the table below lists; one numerical example of the paper is replaced (finding E6).

Every result of the paper is proved for Gaussian noise: the rank-one results of Sections 3
to 5, the estimator of Section 8, and the rank-`r` results of Section 7 and Appendices D
and E. Five rank-one results are also proved at the paper's own noise class, any fixed
entry law with mean 0, variance 1 and a finite fourth moment: the single-table law,
delocalization, unweighted stackSVD, unweighted SVDstack and the estimator. These 5 need
one condition the paper does not state, a density of the entry law. Two remarks of the
paper have no Lean form; "What is not proved" lists them with the other caveats.

The formalization produced 12 findings about the paper. They are hypotheses that a
statement needs (`M ≥ 2`, some `θ_i ≠ 0`, an ordering of the spikes), a numerical example
that sits exactly on its threshold, a proof branch with no argument, 2 results that hold in
a stronger form, and errata. They are proposals to the authors and none is applied to the
paper; `docs/THEOREMS.md` marks each one (E1 to E12) where it arises.

## Start here

1. A reader without Lean: `GETTING_STARTED.md` (by Phillip Nicol) walks through the install,
   the build and one `#check`. A fresh clone builds in 13 minutes on 6 cores and takes 10 GB
   of disk.
2. A coauthor or a statistician: this page, then `docs/THEOREMS.md`, which puts the paper's
   LaTeX statement and the Lean statement side by side for every result and explains each
   hypothesis in words.
3. An auditor: `docs/AUDIT_DOC.md`, which is self-contained: what is claimed, what is
   assumed, how each claim is checked, and the transcript of the build and of the 8 gate
   scripts.

## What is proved

The model of every row is the paper's: `M` tables `X_i = θ_i u_i vᵀ + Z_i / √d`, one shared
direction `v`, entries of `Z_i` i.i.d. `N(0, 1)`, tables independent, `M` fixed while
`n_i, d → ∞` with `n_i / d → c_i > 0`. "Performance" is `⟨v̂, v⟩²` for the estimate `v̂`.
Every limit is convergence in probability, as in the paper. `β_i² = (θ_i⁴ − c_i)/(θ_i⁴ + θ_i²)`
when `θ_i⁴ > c_i` and `0` otherwise. "Extra assumptions" lists what a row needs beyond this
model; "none" means the paper's own hypotheses only. Each Lean name is a theorem in
`lean/StackedSVD/<file>`; the `_gaussian` suffix marks the version whose only inputs are the
model hypotheses, and the `_of_moments` suffix the version at four moments plus a density
(caveat 1).

| Paper result | What is proved, in words | Extra assumptions | Lean theorem (file) |
|---|---|---|---|
| `prop:single_table` | One table: `⟨v̂, v⟩² → β²`, and the top eigenvalue of `XᵀX` converges to `θ² + 1 + c + c/θ²` above the threshold `θ⁴ > c` and to the bulk edge `(1 + √c)²` below it | none | `singleTableLaw_of_gaussian` (`RMT/Full.lean`); `singleTableLaw_of_moments` (`RMT/General/Edge/Sup.lean`) |
| `lem:delocalization` | Per-table estimates of two different tables, each signed so that `⟨v̂_i, v⟩ ≥ 0`: `⟨v̂_i, v̂_j⟩ → β_i β_j` | none | `lem_delocalization` (`SVDStack/Gram.lean`) with the law discharged by `singleTableLaw_of_gaussian`; `lem_delocalization_of_moments` (`General/Layer1.lean`) |
| `prop:stacksvd_general` | Unweighted stackSVD: `⟨v̂, v⟩² → ((Σθ_i²)² − Σc_i) / ((Σθ_i²)(Σθ_i² + 1))` above the threshold `(Σθ_i²)² > Σc_i`, `0` below | none | `prop_stacksvd_general_gaussian` (`StackSVD/Main.lean`); `prop_stacksvd_general_of_moments` (`General/Layer1.lean`) |
| `thm:svd_stack_general` | Unweighted SVDstack: `⟨v̂, v⟩² → (βᵀ v_max(A_β))² / λ_max(A_β)` with `A_β = ββᵀ + diag(1 − β_i²)` | at least 2 tables with `β_i > 0` (caveat 3) | `thm_svd_stack_general_gaussian` (`SVDStack/Main.lean`); `thm_svd_stack_general_of_moments` (`General/Layer1.lean`) |
| `thm:simple_thm1` (cor. 1) | Identical tables (`θ_i = θ₀`, `c_i = c₀`): both closed forms, both branches of the threshold, every `M ≥ 1` | none | `thm_simple_thm1_stacksvd_gaussian` (`StackSVD/Main.lean`), `thm_simple_thm1_svdstack_gaussian_full` (`SVDStack/Simple.lean`) |
| `cor.2` (binary weights) | Keep any nonempty subset `S` of tables: the limit is the stackSVD formula on `S`; the best subset exists | none | `stackPerfW_binary_tendsto_gaussian`, `exists_binary_tendsto_max_gaussian` (`StackSVD/Weighted.lean`) |
| `thm:stacksvd_weighted` | Weights `w_i ∝ θ_i / √(θ_i² + c_i)` are optimal; the limit `γ⋆` solves `Σ θ_i⁴ (1 − x)/(c_i + xθ_i²) = 1` when `Σ θ_i⁴/c_i > 1`, else `0`; every other weighting does no better | some `θ_i ≠ 0` (caveat 4) | `thm_stacksvd_weighted_gaussian`, `_opt` (`RMT/Het/Sup.lean`) |
| `thm:svdstack_weighted` | Weights `w_i = θ_i √((θ_i² + 1)/(θ_i² + c_i)) 1{θ_i⁴ > c_i}` give `⟨v̂, v⟩² → S/(S + 1)`, `S = Σ β_i²/(1 − β_i²)`; no nonzero weighting does better (stronger than the paper, finding E1) | at least 1 table with `β_i > 0` | `thm_svdstack_weighted_gaussian`, `_opt_full` (`SVDStack/Weighted.lean`, `SVDStack/Rayleigh.lean`) |
| `lem:secular_equation` | The determinant identity for a diagonal-plus-rank-one matrix, and: a root of the secular equation above the diagonal is the largest eigenvalue, and that eigenvalue is simple (finite matrix algebra, no limit) | none | `det_Rmat_sub`, `lamMax_Rmat_eq`, `topSimple_Rmat` (`Secular.lean`) |
| `thm:stacksvd_binary_optimal_svd_stack` | Binary stackSVD (keep the detectable tables) is at least as good as optimally weighted SVDstack; strictly better with 2 detectable tables | `c_i ≤ 1` on the kept tables only (weaker than the paper) | `thm_stacksvd_binary_optimal_svd_stack_gaussian`, `_strict` (`StackSVD/Weighted.lean`, `StrictFacades.lean`) |
| `prop:dominance` | Optimally weighted stackSVD is at least as good as unweighted stackSVD and as optimally weighted SVDstack; strict under the paper's two side conditions | some `θ_i ≠ 0` (caveat 4) | `prop_dominance_gaussian`, `_strict` (`RMT/Het/Sup.lean`, `StrictFacades.lean`) |
| `prop:binarystacksvd_inadmissable` | For every `ε ∈ (0, 1)` there is an instance with `M = ⌈e^{−γ} e^{2/ε}⌉` tables where optimally weighted stackSVD exceeds `1 − ε` while every binary weighting, optimally weighted SVDstack and both unweighted methods have limit `0`; the instance is built, not assumed | none | `prop_binarystacksvd_inadmissable_exists` (`Existence.lean`), the clauses in `RMT/Het/Sup.lean`, `SVDStack/Inad.lean` |
| `thm:theta_est` | Two tables, `θ_1⁴ > c_1`: the estimator `θ̂_2` of `eq:theta_estimation` is consistent | none | `thm_theta_est_gaussian` (`ThetaEst.lean`); `thm_theta_est_of_moments` (`General/Layer1.lean`) |
| `app:wstacksvd_mle` | The weighted stackSVD objective is the marginal Gaussian log-likelihood of the random-effects model, at finite `N`; its maximizers are the top eigenvectors | none (deterministic identity) | `mleLogLik_eq`, `mleLogLik_max_iff_mem_topSpace`, `thm_wstacksvd_mle_marginal` (`MLE.lean`, `MLEConverse.lean`, `MLEMarginal/Main.lean`) |
| `remark:stack_outperform_svd` | An explicit family (every `M ≥ 2`, all `θ_i = c_i = 1`) where unweighted stackSVD has limit `1 − 2/(M + 1)` and SVDstack has limit `0` under every weighting; the model is built | none | `remark_stack_outperform_svd_exists` (`Existence.lean`, `Remarks.lean`) |
| `remark:svd_outperform_stack` | An explicit 2-table instance (`θ = (√5, 4)`, `c = (1, 38.4)`) where unweighted SVDstack has limit `8/9` and every binary stackSVD at most `2008/2310`; the paper's 3-table example sits on the threshold (finding E6) | none | `remark_svd_outperform_stack_exists` (`Existence.lean`, `Remarks.lean`) |
| Section 7 and Appendix E (rank `r`) | `lem:general_rank_delocalization`, `prop:general_rank_unweighted_svdstack`, `prop:stacksvd_subspace`, `thm:gen_rank_weight_svdstak` (uniform in the weights, finding E2), the worked example, `thm:rank_r_svdstack` (both clauses), `thm:rank_r_stacksvd` | spikes inside a table strictly decreasing and nonnegative (caveat 5) | `*_gaussian` theorems in `RankR/GeneralGaussian.lean`, `RankR/SubspaceGaussian.lean`, `RankR/WeightedUpperG.lean`, `RankR/AlignedMain.lean`, `RankR/Het/Sup.lean` |
| `prop:gen_rank_stacksvd_singleweight` (Appendix D) | StackSVD with one weight `w_i` per table, general rotations `R_i`: the limit of `‖V̂ᵀ V‖_F²` is the paper's sum over the `r` roots `γ_ℓ` of a matrix secular equation, and per pair `(ℓ, k)` the overlap `⟨v̂_ℓ, v_k⟩²` has the limit `swTerm_ℓ (z_ℓ)_k²` | the paper's separation assumption (distinct roots above the threshold, as data `γ`, `z`); no condition on the weights | `prop_gen_rank_stacksvd_singleweight_gaussian`, `_inner_gaussian` (`RankR/SingleWeight/Het/Sup.lean`) |
| `prop:singleweight_suboptimality` (Appendix D) | A built 2-table rank-2 instance (`θ = 8/5`, `c = 1`, one table rotated) where unweighted SVDstack has limit `2β²` and stackSVD with one weight per table has a smaller limit at every positive weight pair | none | `prop_singleweight_suboptimality_gaussian` (`RankR/SingleWeight/Suboptimality.lean`) |

The paper cites 2 limit laws (Liu et al. 2023; Ding 2020) as black boxes. The Lean
development proves them for Gaussian noise from scratch. Each result exists in 2 forms: a
"Layer 1" theorem that takes the limit law as an explicit hypothesis, and the `_gaussian`
theorem above, where that hypothesis is discharged. The Gaussian core is self-contained (95
modules, 75 of them ours, 42494 lines) and imports nothing about stacking or about the
paper, so it can be reused alone: import `StackedSVD.RMT.Full` and read
`SpikedModel.singleTableLaw_of_gaussian`. `docs/TECHNICAL.md` describes the core;
`docs/THEOREMS.md` section 2 lists every hypothesis structure and the theorem that
discharges it.

## What is not proved

1. Noise. Every `_gaussian` theorem assumes Gaussian noise; the paper allows any noise with
   unit second moment and a finite fourth moment. The 5 `_of_moments` theorems of the table
   hold at that class (`NoiseLaw`, `Prob/NoiseLaw.lean`), plus a Lebesgue density of the
   entry law, which the paper does not assume. A law without a density (Rademacher) is
   open. The weighted and the rank-`r` limit laws are proved for Gaussian noise only, so
   every other row is Gaussian-only.
2. Rates. Every limit is convergence in probability, as in the paper. No rates.
3. One detectable table. Unweighted SVDstack (`thm:svd_stack_general`) needs 2 tables with
   `β_i > 0`. With exactly 1, `A_β` is the identity and the paper's formula is not defined;
   the paper leaves that case open.
4. `θ ≡ 0`. `thm:stacksvd_weighted` and `prop:dominance` need some `θ_i ≠ 0`. At `θ ≡ 0`
   the stated optimal weights are all zero and the claim is false (finding E8).
5. Rank `r`. The spikes inside a table are strictly decreasing and nonnegative, a model
   field, so `ℓ_j = j` in `thm:rank_r_stacksvd`; the paper's index `ℓ_j` at unordered `θ` is
   finding E7. That theorem also takes `hnz` (every component has a positive spike in some
   table), which the paper's own hypothesis implies for `r ≥ 2`.
6. Two remarks. The removal remark after `thm:gen_rank_weight_svdstak` (paper line 901) is
   wrong in general and has no Lean form (finding E2). The `M → ∞` remark (line 407) has
   none either: `M` is fixed in every statement.

## Check it yourself

Each Lean name of the table opens in a text editor: `lean/StackedSVD/<file>`, search the
name, read the lines between `theorem <name>` and `:=`. The hypotheses are the binders; the
conclusion follows the last `:`. `lean/StackedSVD/Main.lean` (40 theorems, one per paper
result) shows every statement at once. A theorem is only as good as the definitions in its
statement. Those sit in 5 files, about 1.7k lines (`Defs.lean`, `StackSVD.lean` with
`StackSVDWeighted.lean`, `SVDStack/Defs.lean`, `RankR/General.lean`); the notation table
of `docs/TECHNICAL.md` maps the paper's symbols to them.

To rerun the checker on your own machine (`GETTING_STARTED.md` explains each step):

```sh
curl https://elan.lean-lang.org/elan-init.sh -sSf | sh   # once: elan, the Lean version manager (or: brew install elan-init)
git clone https://github.com/phillipnicol/stackedSVD.git
cd stackedSVD/lean_formalization/lean
lake exe cache get    # once: downloads the prebuilt Mathlib (about 1 GB, 7 GB unpacked)
lake build            # rechecks every proof; must end with "Build completed successfully"
cd ..
scripts/check_sorries.sh                                  # gate 1: no open proof
scripts/check_axioms.sh                                   # gate 2: only the three standard axioms
python3 scripts/check_root_imports.py                     # gate 3: every module is in the audited set
python3 scripts/check_theorems_sigs.py docs/THEOREMS.md   # gate 4: the quoted statements match the source
scripts/check_layering.sh                                 # gate 5: no Gaussian result uses a limit-law hypothesis
python3 scripts/check_paper_edits.py                      # gate 6: prints SKIP here (needs the paper snapshot; the result of record is in docs/TECHNICAL.md)
scripts/check_kernel.sh                                   # gate 7: every declaration re-checked by the kernel (12 min, 20 GB of RAM with 3 workers; KERNEL_CHECK_JOBS=2 for less)
python3 scripts/check_core_imports.py                     # gate 8: the Gaussian RMT core imports no paper-specific module
```

`lake build` is the core check: Lean recompiles and rechecks every proof, and an open proof
would print `declaration uses 'sorry'`. `scripts/README.md` explains what each gate proves.
The gates ran on Linux, by hand; there is no continuous integration.

## Build of record

Lean tree of commit `416ab9d` (2026-09-10): `Build completed successfully (8978 jobs)`; 0
tracked `sorry` sites; 5631 declarations audited, no axiom outside `propext`,
`Classical.choice` and `Quot.sound`; 148 signatures quoted in `docs/THEOREMS.md` match the
source; 231 modules replayed through the kernel, 0 problems; our 210 files compile with 0
linter warnings. `docs/TECHNICAL.md` holds the build of record with its timings, and
`docs/AUDIT_DOC.md` Appendix A the full transcript of the 8 gates.

## Layout

| Path | What it is |
|---|---|
| `README.md` | this page |
| `GETTING_STARTED.md` | install Lean, build the project, check one theorem yourself |
| `docs/` | the detailed documents; `docs/README.md` names each one (`THEOREMS.md`, `AUDIT_DOC.md`, `TECHNICAL.md`, `SORRIES.md`) |
| `lean/` | the Lean project (Lean 4 v4.33.0, Mathlib v4.33.0 and StatsMLlib `37286c3`, pinned in `lakefile.toml` and `lean-toolchain`); the source is `lean/StackedSVD/` (108508 lines, 210 files, plus 20 vendored files); `lean/StackedSVD/Main.lean` states the 40 theorems |
| `scripts/` | the 8 gate scripts and the numeric checks; `scripts/README.md` describes each |
| `notes/` | the lab notebook of the development repository; not part of this copy (see "About this copy" below) |
| `CITATION.cff` | how to cite the paper and this repository |
| `LICENSE` | Apache License 2.0 |

## How it was made

Claude (Anthropic) wrote the Lean statements and proofs between 2026-08-29 and 2026-09-10,
under a fixed workflow: for each result the authors reviewed the Lean statement against the
paper before a proof started, and 8 mechanical gates ran before each push. While this work
was under way, Anthropic published "Formalizing Fermat's Last Theorem"
(<https://www.anthropic.com/research/formalizing-fermats-last-theorem>, 2026-09-04), a
report that Claude wrote a complete Lean proof of Fermat's Last Theorem in 11 days. On
2026-09-08 OpenAI published "On the Navier-Stokes Millennium Prize Problem"
(<https://openai.com/index/navier-stokes-solution/>), an AI-generated proof that the forced
three-dimensional equations develop a singularity in finite time, with a Lean
formalization. Each of those works formalizes one theorem. The paper here states many
separate limits and comparisons, and the work went as much into stating each one
faithfully as into proving it.

## About this copy

This folder is the public copy of the development repository of the formalization, exported
on 2026-09-11 from commit `41e730f` by its script `scripts/release/make_release.py`. The
Lean sources, the gate scripts and the numeric scripts are byte-identical to that commit. Not
part of this copy: `notes/` (the lab notebook of the development, 374 files: the review
note of every result, the audit packets, the campaign logs, the read-only snapshot of the
paper, and the two paper-facing documents, private until the revised paper is public),
`CLAUDE.md` (the working rules of the AI assistant) and the scratch folder
of the count proof. The documents and the Lean docstrings cite those files by path;
`scripts/release_withheld.txt` lists every withheld path (449 in all), and
`scripts/check_paths.py` reads it, so a citation of a withheld file counts as withheld, not as
missing. Gate 6 (`check_paper_edits.py`) needs the paper and prints `SKIP` here; its result of
record is in `docs/AUDIT_DOC.md`, Appendix A.

## Questions, citation and license

Questions and corrections: open an issue on the GitHub repository, or write to the authors
(the addresses are in the paper). `CITATION.cff` gives the citation of the paper
(arXiv:2507.22170) and of this repository. The license is Apache 2.0 (`LICENSE`), the license
of Lean, Mathlib and StatsMLlib. The 20 files under `lean/StackedSVD/Vendor/COLT83/` are
backported Mathlib work by Rémy Degenne (copyright the author, the same license); their
headers and the `-- Backported` comments mark the origin and the changes.
