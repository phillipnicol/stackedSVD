# TECHNICAL.md: the detailed map of the formalization

> Location (2026-09-03): this file lives in `docs/`. The bare names `TECHNICAL.md`, `THEOREMS.md`
> and `AUDIT_DOC.md` refer to the files of this folder; every other path is relative to the
> repository root. `paper_edits.md` and `NOTE_FOR_COAUTHORS.md` moved from this folder to
> `notes/` on 2026-09-11; this file names them `notes/paper_edits.md` and
> `notes/NOTE_FOR_COAUTHORS.md` throughout.
> The short `README.md` at the root is the plain-language entry point.

Lean 4 formalization of *"Stacked SVD or SVD stacked? A random matrix theory perspective on
data integration"* (Baharav et al., AoS resubmission; arXiv:2507.22170 is the older public
version with the same labels). It builds on Mathlib **v4.33.0** and StatsMLlib `37286c3`,
plus a small set of vendored Gaussian-comparison files (COLT83, Apache-2.0). The tree has
**no open proof** since 2026-09-10 evening, when the last tracked count, the Furedi-Komlos
count now in `RMT/General/Edge/CountBound.lean` (`cellCard_mul_le_pos`), was proved, and it
contains no `axiom`. Every public declaration depends only on `propext`, `Classical.choice`
and `Quot.sound`; no declaration depends on `sorryAx`. Gate 2 (`scripts/check_axioms.sh`)
allows `sorryAx` only for a declaration whose `sorry` has a row in `docs/SORRIES.md`, and
that table is empty.
Build of record: root `lake build` of 2026-09-10 15:58:20 to 16:01:47 EDT
(`lake-build-capped.sh 8`) on the working tree of commit `416ab9d` on `stage3` (the tree is
exactly that commit's tree; nothing was uncommitted except the scratch folders, which are
not in the build): **8978 jobs, exit 0**, 0 errors, 6 jobs built (`RMT.General.Edge.CountBound`,
`Edge.Excess`, `Edge.Markov`, `Edge.Sup`, `General.Layer1`, the root `StackedSVD`), the rest
replayed (built by the root build of 15:53:55 to 15:57:59 EDT on the same sources up to one
docstring line, exit 0, 8978 jobs, which compiled `Edge.Dyck` in 42 s, `Edge.Code` in 58 s
and `Edge.CountBound` in 32 s, and by the integration unit's builds of 15:40 to 15:46 EDT);
0 linter warnings in our 210 files, the 20 vendored files keep their 10
header-linter lines. Eight gates on that tree, each exit 0 (the seven quick ones at 16:02
EDT, gate 7 at 16:04:23 to 16:15:50 EDT): `check_sorries.sh` 0 sorry sites in 0
declarations, 0 rows, all matched; `check_axioms.sh` **5631 declarations** audited (5631
3297: the auditor's `AXIOMAUDIT_END <reported> <skipped>` counts, 5631 declarations reported
and 3297 compiler-generated auxiliaries skipped, `scripts/AxiomAudit.lean`), no axiom outside `propext`, `Classical.choice` and `Quot.sound` and no untracked
sorry; no declaration depends on `sorryAx`; `check_root_imports.py`
230 modules, 210 ours, 229 reachable; `check_theorems_sigs.py` 148
signatures, 136 verbatim, 12 documented abbreviations; `check_layering.sh` 3 endpoints, 0
forbidden hits; `check_paper_edits.py` 29 replacements, 6 amended blocks;
`check_core_imports.py` (v2) 95 core modules (75 ours, 42494 lines), no import outside the
core, 177 EXT modules (157 ours, 87971 lines), import-closed, 5 entry theorems present;
`check_kernel.sh` 231 modules replayed through the kernel, 0 problems, every one of the 231
oleans of the package (686 s, 3 workers; rerun on the same tree at 18:14:24 to 18:25:55 EDT
after the script gained its portable mode, 690 s, the same verdict).
`scripts/check_paths.py` (not a gate of record): 2087 paths resolve, 0 missing (2026-09-10 evening, after the release edits). Earlier
builds of record (every one since 2026-09-02, the text unchanged): `notes/BUILD_HISTORY.md`.

`THEOREMS.md` is the self-contained auditor document: paper statement, Lean statement,
hypotheses in words, modeling choices, scope notes, and audit trail for every result below.
`notes/paper_edits.md` carries the twelve findings E1 to E12 that the formalization produced about the
paper itself.

## Notation, paper to Lean

Every Lean name below lives in `lean/StackedSVD/`; `grep -n "def <name>"` or
`grep -n "structure <name>"` in the file finds it (line numbers are not quoted here because
they drift). `m` is the model, `N` the index of the sequence, `ω` the sample point. This
table was the "Notation" section of `README.md` until 2026-09-10.

| Paper | Lean | File |
|---|---|---|
| `M` tables of size `n_i × d`, `n_i / d → c_i` | `M : ℕ`, `n : Fin M → ℕ → ℕ`, `d : ℕ → ℕ`, `(m.tbl i).Regime (c i)` | `Defs.lean` (`SpikedModel.Regime`) |
| one table `X = θ u vᵀ + Z / √d` | `SpikedModel` (fields `θ`, `u`, `v`, `Z`), the data `m.X N ω`, the scaled noise `m.E N ω` | `Defs.lean` (`SpikedModel`, `SpikedModel.X`, `SpikedModel.E`) |
| `M` tables with one shared `v` | `MultiTableModel` (`tbl`, `hv`) | `Defs.lean` |
| `Z_i` with i.i.d. `N(0, 1)` entries, tables independent | `gaussianMatrix`, `SpikedModel.GaussianNoise`, `MultiTableModel.JointGaussianNoise` | `Defs.lean` |
| performance `⟨v̂, v⟩²` of the top right singular vector of `X` | `overlap X v` (`= ‖P v‖²` for the top spectral projector `P` of `Xᵀ X`) | `Defs.lean` |
| `β_i²`, `ρ²`, the bulk edge `(1 + √c)²` | `betaSq θ c`, `rhoSq θ c`, `bulkEdge c` | `Defs.lean` |
| convergence in probability | `TendstoInProb μ f a` | `Defs.lean` |
| the stack `[X_1; …; X_M]`, its limit `((Σθ_i²)² − Σc_i) / ((Σθ_i²)(Σθ_i² + 1))` | `m.stack` (one `SpikedModel` with `θ = ‖θ‖₂`), `stackSVDLimit θ c` | `StackSVD.lean` |
| the weighted stack `[w_1 X_1; …; w_M X_M]` and its performance | `m.stackW w`, `m.stackPerfW w N ω` | `StackSVDWeighted.lean` |
| `A_β`, `S = Σ β_i² / (1 − β_i²)`, `w⋆`, the svdstack limits | `Abeta β`, `Sval β`, `optW β`, `svdstackLimit β`, `svdstackLimitOpt β` | `SVDStack/Defs.lean` |
| svdstack: the top right singular vector of `[v̂_1ᵀ; …; v̂_Mᵀ]`, and its performance | `m.svdstackEst N ω`, `m.svdstackPerf N ω` | `SVDStack/Defs.lean` |
| rank `r_i` table `X_i = U_i Θ_i V_iᵀ + Z_i / √d` | `SpikedModelR` (fields `θ`, `U`, `V`, `Z`) | `RankR/General.lean` |
| shared subspace `V` and alignments `R_i`, `V_i = V R_i` | `UnalignedModelR` (fields `V`, `R`, `hv`) | `RankR/General.lean` |
| rank-`r` performance `‖V̂ᵀ V‖_F²`, unweighted and weighted | `m.perfRG N ω`, `m.perfRGW W N ω` | `RankR/General.lean` |

## Main results

`*_gaussian` names take model hypotheses only (independent Gaussian tables, proportional
regime). The Layer 1 name takes the random matrix theory input as an explicit hypothesis
structure. Every convergence is in probability.

Reading the table. A name that starts with `_` continues the name before it: in the row for
`thm:svd_stack_general` the Layer 1 column reads `thm_svd_stack_general`, `_inner`, `_zero`,
that is three declarations `thm_svd_stack_general`, `thm_svd_stack_general_inner` and
`thm_svd_stack_general_zero`. The file column names the file of the main declaration of the
row; a scalar half in `Scalars.lean` and a Gaussian discharge in `RMT/` or `RankR/RMT/` may
live elsewhere, and `THEOREMS.md` gives the file of every declaration.

A name that ends in `_general` or `_of_general` is a general-noise theorem and sits in the
Gaussian column next to its Gaussian twin. It takes the paper's own noise class
`assum:general_noise`, that is independent tables with i.i.d. entries of one fixed law with
mean 0, variance 1 and a finite fourth moment. It adds two binders the paper does not have,
a Lebesgue density on the law and the upper edge of the bulk. It is therefore an
implication, not an unconditional theorem.
Five rows carry one (Stage 0 of `notes/NONGAUSSIAN_SCOPE.md`, 2026-09-08, and Stage 1,
2026-09-09); `THEOREMS.md` sections 2.1 and 3.1 give both binders and the route.

| Paper label | Gaussian theorem | Layer 1 theorem | File | Scope note |
|---|---|---|---|---|
| `prop:single_table` | `singleTableLaw_of_gaussian`; `singleTableLaw_of_general` (any fixed noise law with mean 0, variance 1, a finite fourth moment and a Lebesgue density, given the upper edge of the bulk, 2026-09-09) | `SpikedModel.SingleTableLaw` (structure) | `RMT/Full.lean`, `RMT.lean`, `RMT/General/Sup.lean` | both regimes; `topSimple` added; the general form is an implication, the edge is a hypothesis |
| `lem:delocalization` | via `lem_delocalization` + `singleTableLaw_of_gaussian` (`RMT/Full.lean`); `lem_delocalization_of_general` (general noise, one edge hypothesis per table, 2026-09-09) | `lem_delocalization` | `SVDStack/Gram.lean`, `General/Layer1.lean` | signed limit `⟪v̂_i, v̂_j⟫ → β_iβ_j` |
| `lem:entrywise_conv_eigenvec` | deterministic | `lem_entrywise_conv_eigenvec`, `lem_entrywise_conv_eigenvec_signed` | `SVDStack/Deterministic.lean`, `SVDStack/EntrywiseSigned.lean` | projector form, sign-free, and the paper's signed coordinate form |
| `prop:stacksvd_general` | `prop_stacksvd_general_gaussian`, `_inner_gaussian`; `prop_stacksvd_general_of_general` (general noise, one edge hypothesis at the stack, 2026-09-09) | `prop_stacksvd_general`, `_inner` | `StackSVD/Main.lean`, `General/Layer1.lean` | reduces to the stack at `(‖θ‖₂, ‖c‖₁)`; `stack_law_general` keeps that reduction at a general law |
| `thm:svd_stack_general` | `thm_svd_stack_general_gaussian`, `thm_svd_stack_general_inner_gaussian`, `thm_svd_stack_general_zero_gaussian`; `thm_svd_stack_general_of_general` (general noise, one edge hypothesis per table, 2026-09-09) | `thm_svd_stack_general`, `_inner`, `_zero` | `SVDStack/Main.lean`, `General/Layer1.lean` | `β_2 > 0` unsorted; the general form covers the main clause, not `_zero` |
| `thm:simple_thm1` (cor. 1) | `thm_simple_thm1_stacksvd_gaussian`, `thm_simple_thm1_svdstack_gaussian_full` | `Scalars.simple_thm1_stacksvd`, `_svdstack` | `StackSVD/Main.lean`, `SVDStack/Simple.lean` | svdstack half for every `M ≥ 1` |
| `cor.2` (binary) | `stackPerfW_binary_tendsto_gaussian`, `exists_binary_tendsto_max_gaussian` | `Scalars.binaryStackSVDLimit_eq_stackSVDLimit` | `StackSVD/Weighted.lean` | any nonempty subset, strict threshold |
| `thm:stacksvd_weighted` | `thm_stacksvd_weighted_gaussian`, `_gaussian_opt`, `_gaussian_inner`, `_gaussian_smul` | `thm_stacksvd_weighted`, `_inner`, `_general` | `RMT/Het/Sup.lean`, `StackSVD/Weighted.lean` | both regimes, exact edge, `∃ i, θ_i ≠ 0` |
| `thm:svdstack_weighted` | `thm_svdstack_weighted_gaussian`, `_gaussian_opt`, `thm_svdstack_weighted_gaussian_opt_full`, `thm_svdstack_weighted_paper_gaussian` | `thm_svdstack_weighted`, `_paper`, `_general`, `_inner`, `_zero` | `SVDStack/Weighted.lean`, `SVDStack/Rayleigh.lean` | `_opt_full` is uniform in `w` (finding E1) |
| `lem:secular_equation` | deterministic | `secular`, `IsGammaTop`, `det_Rmat_sub`, `lamMax_Rmat_eq`, `topSimple_Rmat` | `StackSVDWeighted.lean`, `Secular.lean` | finite matrix algebra |
| `thm:stacksvd_binary_optimal_svd_stack` | `thm_stacksvd_binary_optimal_svd_stack_gaussian`, `thm_stacksvd_binary_optimal_svd_stack_gaussian_strict` | `Scalars.svdstackOpt_le_binary`, `Scalars.svdstackOpt_lt_binary` | `StackSVD/Weighted.lean`, `StrictFacades.lean` | `c_i ≤ 1` on kept tables only |
| `prop:dominance` | `prop_dominance_gaussian`, `prop_dominance_gaussian_strict`, `_strict_svdstack`, `_strict_unweighted` | `Scalars.svdstackOpt_le_stackSVDLimitW`, `stackSVDLimit_le_stackSVDLimitW` | `RMT/Het/Sup.lean`, `StrictFacades.lean` | both dominance clauses |
| `prop:binarystacksvd_inadmissable` | `prop_binarystacksvd_inadmissable_exists` (the paper's existential form, model built), `prop_binarystacksvd_inadmissable_gaussian`, `stackPerfW_binary_inad_gaussian`, `stackPerfW_opt_inad_gaussian`, `svdstackPerf_inad_gaussian`, `svdstackPerfW_opt_inad_gaussian` | `Scalars.binarystacksvd_inadmissable`, `_ceil` | `Existence.lean`, `RMT/Het/Sup.lean`, `SVDStack/Inad.lean`, `Scalars.lean` | explicit instance, uniform over subsets; all five clauses on the model, and the `∃` form at the paper's ceiling |
| `thm:theta_est` | `thm_theta_est_gaussian`; `thm_theta_est_general` (any fixed noise law with mean 0, variance 1, finite fourth moment; the single-table law of table `i` stays a hypothesis, 2026-09-08); `thm_theta_est_of_general` (that law discharged too, at the price of a density and the edge of table `i`, 2026-09-09) | `thm_theta_est` | `ThetaEst.lean`, `General/Layer1.lean` | item P proved, not assumed; `ThetaEstLaw` also discharged at `NoiseLaw` (`thetaEstLaw_of_general`) |
| `app:wstacksvd_mle` | deterministic | `mleLogLik_eq`, `mleLogLik_max_iff_mem_topSpace`, `thm_wstacksvd_mle_marginal` | `MLE.lean`, `MLEConverse.lean`, `MLEMarginal/Main.lean` | finite-`N` identity, its converse, and the Gaussian marginalization from the random-effects model (L5, 2026-09-03) |
| `remark:stack_outperform_svd` | `remark_stack_outperform_svd_exists` (model built, every `M ≥ 2`), `remark_stack_outperform_svd_stack`, `_svdstack`, `_svdstack_uniform` | instances of the closed forms | `Existence.lean`, `Remarks.lean`, `RemarksUniform.lean` | explicit rationals |
| `remark:svd_outperform_stack` | `remark_svd_outperform_stack_exists` (model built, `θ = (√5, 4)`, `c = (1, 38.4)`), `remark_svd_outperform_stack_two`, `_pair_*`, `_three_*` | instances of the closed forms | `Existence.lean`, `Remarks.lean`, `RemarksUniform.lean` | the `M = 3` example is on the threshold (finding E6) |
| `lem:general_rank_delocalization` | via `tableLawR_of_gaussian_rk` (`RankR/RMT/TableLawGaussian.lean`) | `lem_general_rank_delocalization_general`, `lem_general_rank_delocalization` | `RankR/GeneralMain.lean`, `RankR/Unweighted.lean` | off-diagonal half; `IndepNoise` suffices |
| `prop:general_rank_unweighted_svdstack` | `prop_general_rank_unweighted_svdstack_general_gaussian`, `prop_general_rank_unweighted_svdstack_gaussian`, `prop_general_rank_unweighted_svdstack_general_frobenius_eig_gaussian` | `prop_general_rank_unweighted_svdstack_general`, `_frobenius`, `_frobenius_eig` | `RankR/GeneralGaussian.lean`, `RankR/GeneralFrob.lean`; the `r_i = 1` twin in `RankR/Unweighted.lean` | general `r_i`; trace form and Frobenius form |
| `prop:stacksvd_subspace` | `prop_stacksvd_subspace_gaussian`, `prop_stacksvd_subspace_general_gaussian`, `prop_stacksvd_subspace_frobenius_eig_gaussian`, `prop_stacksvd_subspace_general_frobenius_eig_gaussian` | `prop_stacksvd_subspace`, `_general`, `_frobenius`, `_frobenius_eig`, `_general_frobenius`, `_general_frobenius_eig` | `RankR/SubspaceGaussian.lean`, `RankR/SubspaceMain.lean`, `RankR/Frobenius.lean` | no eigengap of `C`; ordered spikes per table; projector-sum form and Frobenius form at both `r_i = 1` and general `r_i` |
| `thm:gen_rank_weight_svdstak` | `thm_gen_rank_weight_svdstak_general_r_paper_gaussian`, `thm_gen_rank_weight_svdstak_general_r_full_gaussian`, `thm_gen_rank_weight_svdstak_gaussian` | `thm_gen_rank_weight_svdstak`, `_opt`, `_max`, `_paper`, `_general_r`, `_full` | `RankR/GeneralGaussian.lean`, `RankR/WeightedMain.lean` | `rank B_R = r` in place of `β_ij > 0`; uniform-in-`W` bound (finding E2) |
| `eq:psi_equation`, `eq:perf_rankr_ex_*` | `example_perfR_tendsto_gaussian`, `example_perfStackR_tendsto_gaussian` | `example_perfR_tendsto` (`RankR/Example.lean`), `example_perfStackR_tendsto` (`RankR/SubspaceMain.lean`) | `RankR/Example.lean`, `RankR/SubspaceGaussian.lean` | the worked example and its kinks |
| `thm:rank_r_svdstack` | `thm_rank_r_svdstack_aggregate_gaussian`, `thm_rank_r_svdstack_component_gaussian`, `_of_model_gaussian` | `thm_rank_r_svdstack_aggregate`, `_component`, `_component_eig` | `RankR/GeneralGaussian.lean`, `RankR/WeightedUpperG.lean`, `RankR/AlignedMain.lean` | both clauses; canonical frame; not a corollary (finding E4) |
| `thm:rank_r_stacksvd` | `thm_rank_r_stacksvd_gaussian`, `thm_rank_r_stacksvd_proj_gaussian`, `_inner_gaussian`, `_frobenius_gaussian` | `thm_rank_r_stacksvd`, `_proj`, `_inner`, `_frobenius` | `RankR/Het/Sup.lean`, `RankR/StackMain.lean` | ordered spikes, `R_i = 1`, so `ℓ_j = j` (finding E7) |
| `prop:gen_rank_stacksvd_singleweight` | `prop_gen_rank_stacksvd_singleweight_gaussian`, `_inner_gaussian`; `singleWeightLaw_of_gaussian` | `prop_gen_rank_stacksvd_singleweight`, `_inner`; `EigSep`, `SingleWeightLaw` | `RankR/SingleWeight/Het/Sup.lean`, `RankR/SingleWeight/Scalars.lean`, `Defs.lean`, `Main.lean` | one weight per table, general `R_i`; `EigSep.sorted` excludes ties (finding E11); no condition on the weights (F18c dropped D37's `w_i ≠ 0` the same day) |
| `prop:singleweight_suboptimality` | `prop_singleweight_suboptimality_gaussian` (F18b, 2026-09-07: `hlaw_witness` covers the tie, the two-root and the one-root weight pairs) | `prop_singleweight_suboptimality_of_law` (`hlaw` the one hypothesis) | `RankR/SingleWeight/Example.lean`, `Existence.lean`, `Suboptimality.lean` | existence instance `θ_0 = 8/5`; every `w_i > 0` stays below `2 β_0²`; exact limit on the one-root pairs where the paper states a bound (E12) |

## How the random matrix theory input enters

**Two layers.** The paper cites `liu2023asymptotic` Theorems 1 and 2 and
`10.3150/19-BEJ1129` Theorem 2.3 as black boxes. No prover reaches those today, so they enter
as `structure`s, one field per fact a proof reads: `SpikedModel.SingleTableLaw`,
`MultiTableModel.HeteroLaw`, `SpikedModelR.TableLawR`, `UnalignedModel.SubspaceLaw` (and
`SubspaceLawG` at general `r_i`), and `UnalignedModelR.HeteroLawR`. Layer 1 derives the
paper's conclusion from the structure, so it holds for any noise law that satisfies it. Each
structure is then discharged for Gaussian noise by a named theorem:
`singleTableLaw_of_gaussian`, `heteroLaw_of_gaussian`, `tableLawR_of_gaussian_rk`,
`subspaceLaw_of_gaussian`, `subspaceLawG_of_gaussian`, `heteroLawR_of_gaussian`. Nothing is
left assumed. Two edge structures, `HeteroEdge` and `HeteroEdgeR`, are likewise proved by a
sharp Sudakov-Fernique bound (`heteroEdge_of_gaussian`, `heteroEdgeR_of_gaussian`).

**The model class.** Noise is Gaussian and independent across tables; the regime is
proportional (`n_i → ∞`, `d → ∞`, `n_i/d → c_i ∈ (0, ∞)`); the asymptotic index is one natural
number `N`. Every conclusion is convergence in probability, never almost sure, which is the
paper's own mode (`main_paper.tex:234` defines `\pto` that way). Overlaps are
stated in projector form, with an inner-product twin named `_inner` for the paper's own
display. Rank-`r` tables carry strictly ordered nonnegative spikes (`SpikedModelR.hθnn`,
`hθanti`; F8 of 2026-09-05 allows a zero last spike).

**The reusable core.** The 95 modules under `singleTableLaw_of_gaussian` and, since Stage 1
(2026-09-09), its general-noise twin `singleTableLaw_of_general` (`Defs`, `Spectral`, `RMT`,
the 19 files of `RMT/` outside `RMT/Het/` and `RMT/General/`, the 22 of `RMT/General/`, the
14 of `RMT/General/Edge/` (Stage 3, 2026-09-10), the 11
of `Prob/`, `LinAlg/Eigen`, `Frame`, `KyFan`, `SpecInvReindex`, `SpecProjPerturb`,
`TopProjPerturb`, and the 20 vendored files of `Vendor/COLT83/`; 75 files and 42494 lines
ours; 57 modules, 37 files and 21088 lines at the creation of the gate on 2026-09-08) import
no module that mentions stacking, weights or the paper's estimators; `scripts/check_core_imports.py` (gate 8,
2026-09-08) keeps it that way, and `README.md` ("The Gaussian RMT core") describes the entry
theorem in words. The other three Gaussian discharge families import the paper's definition
files for the structures they discharge (`RMT/Het/` imports `StackSVDWeighted.lean` for
`HeteroLaw` and `stackW`; `RankR/RMT/`, `RankR/Het/` and `RankR/SingleWeight/Het/` import
`RankR/General.lean`, `RankR/Subspace.lean`, `RankR/StackGamma.lean`, `RankR/StackMain.lean`,
`RankR/SingleWeight/Scalars.lean` and `Main.lean`; `LinAlg/SpecIdx*.lean` and `SpecWindow.lean`
import `SVDStack/Defs.lean` for `vMax`). `python3 scripts/check_core_imports.py --boundary`
prints that list. F28 (done 2026-09-08) moved the Layer 1 theorems that used to sit inside
those definition files into six new files above them (`StackSVD/Main.lean`,
`StackSVD/Weighted.lean`, `RankR/Unweighted.lean`, `RankR/WeightedMain.lean`,
`RankR/SubspaceMain.lean`; `SVDStack/DelocDir.lean` carries consumer-free content moved out
of `SVDStack/Gram.lean`), so a discharge family's remaining imports carry no Layer 1
consumer. This does not put the three families inside the core itself (they still import
the paper's model files); it puts them inside a wider, separately gated set, EXT: the
95-module core, the four families minus their three `Sup` files
(`RMT/Het/Sup.lean`, `RankR/Het/Sup.lean`, `RankR/SingleWeight/Het/Sup.lean`, each a Layer 1
corollary and so a consumer by design), plus 21 definition and lemma modules the families
need. EXT holds 177 modules (139 at its creation on 2026-09-08; the 24 files of Stage 1 and the 14
of Stage 3 joined the core), is import-closed, and no declaration inside it takes one of
the eight law structures (`SingleTableLaw`, `HeteroLaw`, `ThetaEstLaw`, `TableLawR`,
`SubspaceLaw`, `SubspaceLawG`, `HeteroLawR`, `SingleWeightLaw`) as a hypothesis.
`scripts/check_core_imports.py` (gate 8) checks the core, EXT, and that five entry
theorems exist. `vMax` did not move: `SVDStack/Defs.lean` has no consumer of any of the
eight structures, so it already belonged to EXT.

**Not covered.** The paper's `assum:general_noise` is wider than Gaussian. The Layer 1
theorems read the noise only through their hypothesis structures and, in the rank-one
SVDstack family, through `MultiTableModel.IndepNoise` (independent tables with any law;
since 2026-09-02, L4; `thm_theta_est` through `MultiTableModel.ThetaEstLaw` since
2026-09-08, F31), so a four-moment discharge of the structures would give every result
unchanged. Stage 1 (2026-09-09) is that discharge for `SingleTableLaw`
(`singleTableLaw_of_general`, `RMT/General/Sup.lean`), and the four representative corollaries
of `General/Layer1.lean` show it plugging in with no Layer 1 statement changed. It is not
unconditional: it carries the upper edge of the bulk as a hypothesis, and it asks for a
Lebesgue density, which excludes Rademacher entries (item F33). Stage 3 (2026-09-10) proves
that edge at four moments (`lamMax_W0_edge_of_general`, `RMT/General/Edge/Sup.lean`, the
upper half of the Bai-Yin theorem, by the moment method with a truncation), so the density
is the one condition beyond the paper that remains. `singleTableLaw_of_moments` and the four
`*_of_moments` facades of `General/Layer1.lean` are the forms with the edge supplied; the
forms with the edge as a binder keep their names. The one combinatorial count inside the
Stage 3 chain, the Furedi-Komlos count, was proved 2026-09-10 evening (now
`cellCard_mul_le_pos` in `RMT/General/Edge/CountBound.lean`); `docs/SORRIES.md` has no row.
The weighted and rank-`r`
structures stay Gaussian, stages 4 and 5 of the same note. One appendix result is not
formalized: the removal
remark at `main_paper.tex:901`, wrong in general (question Q8). One more appendix result,
`prop:singleweight_suboptimality`, was formalized at Layer 1 only until 2026-09-07, under the
hypothesis `hlaw`, the single-weight limit on its two-table witness at every positive weight
pair; the Gaussian discharge of `SingleWeightLaw` (Track G, 2026-09-05) needs the separation
assumption `EigSep`, which fails on that witness at extreme weight ratios. Item F18b (2026-09-07)
proves `hlaw` on the witness regime by regime (`RankR/SingleWeight/Suboptimality.lean`), so
the proposition is now unconditional for Gaussian noise. In the
rank-`r` results the spike order inside a table is a model field, so a table with a repeated
`θ_ij` is outside the class where the paper allows it. `prop:maxrow` and `cor:domination` are
commented out in the paper source, so nothing formalizes them. See `THEOREMS.md` section 7 for
the full list, including the well-formedness hypotheses that remain.

## Verification

```sh
cd lean && lake exe cache get && lake build   # rechecks every proof
cd .. && scripts/check_sorries.sh                   # tree matches the docs/SORRIES.md ledger
CHECK_AXIOMS_BUILD=0 scripts/check_axioms.sh        # every declaration: 3 standard axioms only
python3 scripts/check_root_imports.py                # every module is inside the audited closure
python3 scripts/check_theorems_sigs.py docs/THEOREMS.md  # every quoted signature matches the tree
scripts/check_layering.sh                           # no Gaussian endpoint uses a Layer 1 theorem
python3 scripts/check_paper_edits.py                # the paper edits are exact; amended statements verbatim
scripts/check_kernel.sh                             # every declaration re-checked by the kernel (leanchecker, about 10 min)
```

(On the development server a capped wrapper, `lake-build-capped.sh 12`, stands in for bare `lake build`.) `scripts/README.md` documents every script.
`scripts/make_audit_pack_header.sh` (server only) runs the build and the first four gates in
one pass and prints a dated transcript; `AUDIT_PACK_JOBS` caps the cores of the build it
starts (default 12).

The kernel gate (`check_kernel.sh`, 2026-09-07) is the second pass of the kernel: the
toolchain's `leanchecker` reads the compiled `.olean` of each of our modules and adds every
constant to the environment through the kernel again, with no elaborator, tactic or `unsafe`
code in between. It catches a declaration that skipped the kernel at build time (a test
module that adds `false : False` with `addDeclCore (doCheck := false)` passes `lake build`
and fails the gate; `scripts/README.md`). It trusts Mathlib and StatsMLlib as compiled, and
it is not an external verifier: the kernel that re-checks is the one of the same toolchain.

The root-imports gate matters for a reason that is easy to miss: `AxiomAudit.lean` audits what
`import StackedSVD` reaches, so a module absent from `lean/StackedSVD.lean` would be
compiled by `lake build` and never audited.

Paste this into a Lean file inside the project to check the headline theorems directly:

```lean
import StackedSVD
#print axioms StackedSVD.SpikedModel.singleTableLaw_of_gaussian
#print axioms StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian
#print axioms StackedSVD.MultiTableModel.thm_svdstack_weighted_gaussian
#print axioms StackedSVD.UnalignedModelR.thm_rank_r_stacksvd_gaussian
```

Each prints `[propext, Classical.choice, Quot.sound]`.

**The hypothesis sets are satisfiable.** `StackedSVD/Sat.lean` builds concrete Gaussian models
and applies the three endpoints above to them with no hypothesis left over
(`Sat.Rank1.sat_stacksvd_weighted`, `Sat.Rank1.sat_svdstack_weighted`,
`Sat.RankR.sat_rank_r_stacksvd`), so none of them is vacuous. The module is inside the audited
closure, so the axiom gate rechecks the witnesses on every run. `THEOREMS.md` section 8.1 gives
the models. `StackedSVD/Existence.lean` (2026-09-03) uses the same constructor to prove the
paper's three existential statements with the model built inside the theorem:
`prop_binarystacksvd_inadmissable_exists` (for any `ε ∈ (0, 1)` there is a problem instance
of size `M = ⌈e^{-γ} e^{2/ε}⌉` with the five limits), `remark_stack_outperform_svd_exists`
and `remark_svd_outperform_stack_exists` (`THEOREMS.md` 5.3 and 5.6).

**The two layers are separated at the declaration level.** `scripts/check_layering.sh`
computes the transitive constant closure of the three Gaussian discharges
(`singleTableLaw_of_gaussian`, `heteroLaw_of_gaussian`, `heteroEdge_of_gaussian`) and fails
if it contains a Layer 1 consumer (a theorem that takes `SingleTableLaw` or `HeteroLaw` as a
hypothesis). Before F28 (2026-09-08) the file layout alone did not show this for the
heteroscedastic chain: three files of `RMT/Het/` imported `StackSVDWeighted.lean`, which held
both `HeteroLaw` and its consumer, `thm_stacksvd_weighted_general`. F28 moved every Layer 1
consumer of a discharge family's definition file into a new file above it
(`StackSVDWeighted.lean`'s consumer is now in `StackSVD/Weighted.lean`), so the file layout
shows the split here too. Three `Sup` files still mix a discharge and a consumer on
purpose, since each proves the paper's final Gaussian theorem right where its hypothesis
structure is discharged; the declaration-level check stays the authority regardless of file
layout (`AUDIT_DOC.md` section 4).
`scripts/AuditSignatures.lean` prints the elaborated signature and the axioms of every
headline theorem. `THEOREMS.md` explains how to read the `check_axioms.sh` output.

## Layout

The tree below names every folder and every Lean file, with what it holds. Paths under
`lean/` are relative to the source root `lean/StackedSVD/`. This public copy ships without
`notes/` (so without the paper snapshot and the two paper-facing documents it holds) and
`CLAUDE.md`; `scripts/release_withheld.txt` lists every withheld path.

```
README.md               the plain-language entry point for a coauthor without Lean: the result
                        table, the caveats, how to check it, the folder map
GETTING_STARTED.md      install Lean, build the project, check one theorem; by Phillip
                        Nicol (2026-09-09, as example.md), at the top level since 2026-09-10
CLAUDE.md, LICENSE (Apache 2.0, user choice 2026-09-07), CITATION.cff (the paper and
                        the repository, 2026-09-10)
docs/                   the detailed documents (moved from the top level 2026-09-03):
  README.md             one entry per file of this folder
  TECHNICAL.md          the detailed map: the notation table (paper symbol to Lean name,
                        moved from README.md 2026-09-10), main-results table, the two
                        layers, verification, layout, build of record (the top-level README
                        until 2026-09-03)
  THEOREMS.md           the self-contained auditor annex: paper statement, Lean statement,
                        hypotheses in words, modeling choices, scope, audit trail
  AUDIT_DOC.md          the self-contained audit entry point, with the gate transcript
  SORRIES.md            the ledger of open proofs; empty since 2026-09-10 evening (the last
                        row, the Furedi-Komlos count, was proved then), and gate 1 keeps the
                        ledger and the tree in step
notes/                  living files only: PROGRESS, FLAGGED, FOLLOWUP_LIST, REVISION_LIST,
                        SERVER_TODO, INTERFACES, README (the inventory of paper results
                        and their formalization class), SERVER_SETUP (the toolchain of the
                        shared server and of the laptop), RANK_R_PLAN,
                        NONGAUSSIAN_SCOPE (the stages of the non-Gaussian extension),
                        stage3_edge (the Bai-Yin edge plan; proved 2026-09-10),
                        STAGE3_CAMPAIGN (the decisions, the unit table and the log of the
                        Stage 3 run of 2026-09-10),
                        REPO_CLEANUP_PLAN (the cleanup units of 2026-09-10),
                        BUILD_HISTORY (every earlier build of record, moved out of this file
                        and docs/TECHNICAL.md on 2026-09-10)
notes/archive/          completed material: per-theorem review notes, audits, packs,
                        external reviews and responses, campaign plans, agent reports,
                        PLAN_2026-08-29.md, lean_refs.txt; INDEX.md lists every file
notes/PROGRESS.md       the log; its "Where we are" block is the fastest catch-up
notes/SERVER_TODO.md    per-module compile status and what comes next
notes/FLAGGED.md        decisions Claude took alone (D1 to D38) and questions Q1 to Q17
notes/FOLLOWUP_LIST.md  follow-up items F1 to F31 with their status
notes/REVISION_LIST.md  the paper revisions R1 to R27 that wait for the user
notes/INTERFACES.md     fixed names shared by every note and every Lean file
notes/RANK_R_PLAN.md    Section 7 plan; the code is in lean/StackedSVD/RankR/
notes/paper/main_paper.tex  read-only snapshot of the paper (2026-08-29); never edited (moved
                        from docs/ 2026-09-11)
notes/paper_edits.md   findings E1 to E12 about the paper, each with a Lean anchor and the
                        exact old and new LaTeX of every proposed edit (gate 6), moved from
                        docs/ 2026-09-11
notes/NOTE_FOR_COAUTHORS.md  one page for a coauthor: proved, broke and changed, left; moved
                        from docs/ 2026-09-11
notes/archive/agent_reports/    one report per proof agent
notes/audit_packets/    the external audit packets (rank1_2026-09-03), generated by
                        scripts/make_audit_packet_rank1.py; README.md says what a packet is
notes/discoveries/      short records of things found on the way; README.md lists them
scripts/                README.md (one entry per script), check_sorries.sh, check_axioms.sh, check_root_imports.py,
                        check_theorems_sigs.py, check_layering.sh, check_paper_edits.py,
                        check_kernel.sh (gate 7, 2026-09-07: leanchecker replays every
                        module through the kernel, about 10 min and 15 GB, on the root-disk
                        olean copy after a cmp against the build; since 2026-09-10 evening,
                        on a machine without the nprocs shim, it runs KernelReplay.lean
                        instead, a copy of the toolchain's LeanChecker.lean with a worker
                        pool bounded by -j: 230 modules in 737 s and 20 GB with 3 workers in
                        a fresh clone of the public copy),
                        check_core_imports.py (gate 8, 2026-09-08, widened the same day for
                        F28: the 57 modules of the Gaussian RMT core import no paper-specific
                        module; a wider set, EXT (139 modules: the core, the four Gaussian
                        discharge families minus their three Sup files, and 21 definition and
                        lemma modules), is import-closed and holds no consumer of the eight
                        law structures; five entry theorems are present; under 1 s): the
                        gates before a commit; scripts/numeric/check_singleweight_regimes.py,
                        the numeric check of the three weight regimes of F18b (one weight
                        per table);
                        check_paths.py (2026-09-10, not a gate of record: every path the
                        reader-facing documents name exists; since 2026-09-10 evening it
                        reads scripts/release_withheld.txt in the public copy);
                        release/make_release.py (2026-09-10 evening: exports the public
                        copy, the lean_formalization/ folder of the paper's code repository,
                        notes/RELEASE_PLAN.md); scripts/numeric/ (the
                        top-level scratchpad/ until 2026-09-09: exploration scripts, with
                        their own README)
lean/                   the lake project (Mathlib v4.33.0 + StatsMLlib 37286c3), renamed
                        from StackedSVD/ on 2026-09-10; the sources are lean/StackedSVD/,
                        108508 lines ours in 210 files, 3716 vendored in 20; every path
                        below is relative to lean/StackedSVD/

  Main.lean             the statement file (2026-09-07): 40 theorems, one per paper result
                        in the order of docs/THEOREMS.md, each the final Gaussian (or
                        deterministic) form, proved by one call to the tree; imported last
                        by the root
  Defs.lean             model v2: SpikedModel, MultiTableModel, gaussianMatrix, betaSq,
                        rhoSq, bulkEdge, topProj, overlap, TendstoInProb
  Spectral.lean         overlap bounds, measurability, reindexing, scaling
  RMT.lean              SingleTableLaw: the one hypothesis structure Layer 1 takes
  StackSVD.lean         stacking lemmas, the stack estimator (its Layer 1 theorems moved to
                        StackSVD/Main.lean, F28, 2026-09-08)
  StackSVDWeighted.lean HeteroLaw, secular scalars, optimal weights (its Layer 1 theorems
                        moved to StackSVD/Weighted.lean, F28)
  StackSVD/Main.lean            prop:stacksvd_general (+ Gaussian corollary), thm:simple_thm1
                                stacksvd half (F28, 2026-09-08)
  StackSVD/Weighted.lean        thm:stacksvd_weighted, cor.2 binary,
                                thm:stacksvd_binary_optimal_svd_stack (F28)
  Scalars.lean          class (A): cor.1, cor.2, binary optimal, prop:dominance
  SVDStack.lean         entry point of SVDStack/
  SVDStack/Defs.lean            A_beta, svdstack estimator, Sval, optW
  SVDStack/Deterministic.lean   abeta_gap, lem:entrywise_conv_eigenvec, limit lemmas
  SVDStack/Gram.lean            align, gram, item D lem:delocalization (the consumer-free
                                half moved to SVDStack/DelocDir.lean, F28)
  SVDStack/DelocDir.lean        perpOf, delocDir, measure_deloc_le_of_pi: the consumer-free
                                half of item D, moved out of Gram.lean (F28, 2026-09-08)
  SVDStack/Main.lean            thm:svd_stack_general (projector, inner, zero, Gaussian)
  SVDStack/Weighted.lean        thm:svdstack_weighted and the paper-literal weights
  LinAlg/TopProjPerturb.lean    Weyl and Davis-Kahan for the top projector
  LinAlg/SpecProjPerturb.lean   the same for the top-r spectral projector (rank-r)
  Prob/TendstoInProb.lean       closure lemmas for convergence in probability
  Prob/PolynomialNull.lean      null zero sets of polynomials; gaussianMatrix << volume
  Prob/GaussianMatrix.lean      flattening, transpose and rotation laws
  Prob/GaussianAdapters.lean    Stein identity, Lipschitz concentration
  Prob/WithDensityPi.lean       densities of product measures and of linear images (L5)
  Prob/GaussianDensity.lean     the density of N(0, S) on Fin d → ℝ, image of the standard
                                Gaussian under a matrix, its charFun (L5)
  Prob/NoiseLaw.lean            NoiseLaw (a fixed law: probability, mean 0, second moment 1,
                                finite fourth moment), noiseMatrix, GeneralNoise,
                                JointGeneralNoise, noiseLaw_gaussian, PiLaw.Centered (Stage 0
                                of notes/NONGAUSSIAN_SCOPE.md, 2026-09-08)
  Prob/LinFormMoments.lean      moments of a linear form of i.i.d. coordinates: the fourth
                                moment at most ν₄ Σ a_l⁴ + 3 (Σ a_l²)²
  Prob/NoiseMoments.lean        the second moments of ‖Z a‖²/D and xᵀ Z a under noiseMatrix
  Prob/Chebyshev.lean           the two Chebyshev copies the general chain uses (Stage 1)
  Prob/Tensorization.lean       Efron-Stein tensorization of the variance over Measure.pi
  RMT/MP.lean, MP7.lean         Marchenko-Pastur transform, real and holomorphic branch
  RMT/R0.lean                   splitting bridge X^T X = W0 + q q^T, block laws
  RMT/R4.lean, R4C.lean         deterministic finite-rank layer, complex resolvent
  RMT/R1.lean, ResolvDeriv.lean, SteinStep.lean   complex MP law for the trace
  RMT/R2.lean                   isotropic resolvent forms
  RMT/R3.lean, R3minus.lean     upper and lower edge of the Wishart block
  RMT/T.lean                    transfer from complex z to the real axis
  RMT/R5.lean                   the analytic core (model-free, consumes ResolventLimits)
  RMT/Symmetry.lean             item Sym: delocUniform by rotation invariance
  RMT/R6.lean                   item R6': subcritical align -> 0
  RMT/Simplicity.lean           item S: top eigenvalue simple a.s.
  RMT/Sup.lean, TailShift.lean  milestone L2a and the index shift that drops hn2
  RMT/Full.lean                 singleTableLaw_of_gaussian, both regimes
  RMT/General/                  Stage 1, the non-Gaussian single table (notes/archive/
                                prop_single_table_general.md, 2026-09-09): Defs (the model
                                definitions of the stage and, since F35, the structure
                                ResolventFormsC), QuadForm,
                                Stability, Companion, MPtilde, Trace (the trace law at four
                                moments), R3minus (the lower edge, no hn2), Simplicity, Iso
                                (the isotropic law, a measure bound uniform in the two
                                directions, tendsto_isoRate), IsoMixed (the mixed form
                                xᵀ G Eᵀ y), FormsBridge, ProbC, Forms, FormsLimits, FormsSup,
                                FormsGeneral (resolventFormsC_of_general'), Deloc (the
                                window bound), DelocLimits, DelocAlign
                                (align_tendstoInProb_of_subcritical_general'), DelocUniform (unit
                                G7, delocUniform_of_lamMax: the window bound at w ⊥ v),
                                DelocUniformSup (delocUniform_of_general', both regimes),
                                Sup (resolventLimits_of_general, delocUniform_of_general,
                                singleTableLaw_of_general, the Stage 1 endpoint)
  RMT/General/Edge/             Stage 3, the sharp upper edge at a general law (route C of
                                notes/stage3_edge.md, the moment method with a truncation;
                                the campaign log is notes/STAGE3_CAMPAIGN.md, 2026-09-10).
                                Fourteen files:
    Defs.lean                   the definitions the stage shares: truncMap, the level
                                truncLevel a d = d^(1/2 - a), momOrder, markovConst,
                                fkFactor, TruncNoiseLaw, walkMult, the walk classes
                                IsPaired and IsExcess, excessSum
    Arith.lean                  unit K, first half: 1 <= k_N and 8 k_N <= d_N for N large,
                                fkFactor at most 1/2 for N large, and the Markov limit at
                                C = markovConst c eps
    Trunc.lean                  unit T: the split Z = Zhat + R + m_T J at the level
                                d^(1/2 - a), truncLaw is a TruncNoiseLaw, the mean matrix
                                has small operator norm
    Sparse.lean                 unit L: with probability tending to 1 the discarded part R
                                has operator norm at most eps sqrt d (no large entry, and
                                no row or column with two of them)
    Trace.lean                  unit E: trace ((Y^T Y)^k) as a sum over the closed walks of
                                length 2k on Fin n x Fin d, its integral as a product of
                                entry moments, and the k = 2 closed form as the check
    Compare.lean                unit C: the three walk classes; the trace at a
                                TruncNoiseLaw is at most the Gaussian trace plus excessSum
    Gaussian.lean               unit G: the Gaussian trace bound at k = momOrder C d, from
                                R3.integral_opNorm_le and R3.measure_opNorm_ge_le
    Count.lean                  unit X, part 1: the walk classes and the vertex count
                                (walkVerts, clsCard, cellCard); the Furedi-Komlos count
                                itself moved to CountBound.lean (536 lines, was 612)
    Code.lean                   X8 plan part U: the six-part step code of a walk (code,
                                code_injective), the bad-step bound 4 D by a credit-debt
                                induction (card_codeP_le), the cell code count
                                cellCard_le_code' (2526 lines)
    Dyck.lean                   X8 plan part L: the contour walk of a Dyck word
                                (dyckWalk), the labeled tree-walk lower bound
                                (dyck_label_card_le), and the Narayana lower bound by a
                                double rotation (narayana_lower, tree_cell_lower); 2063
                                lines
    CountBound.lean             X8 plan part A and the assembly: the three arithmetic
                                lemmas (descFac_mul_min_pow_le, choose_mul_le_pow_choose,
                                code_pow_arith), cellCard_mul_le_pos (the Furedi-Komlos
                                count itself, proved 2026-09-10), and cellCard_mul_le,
                                clsCard_mul_le moved here from Count.lean (294 lines)
    Excess.lean                 unit X: excessSum at most 2 fkFactor times the Gaussian
                                trace, relative to the tree-walk class
    Markov.lean                 unit K, second half, and unit M: the two spectral
                                contractions, the rank-one downdate, the bundled moment
                                bound, and the Markov step at the moment order
    Sup.lean                    the assembly: tendsto_measure_opNorm_edge (the canonical
                                law), opNorm_sq_edge_of_general (the target),
                                lamMax_W0_edge_of_general (the hedge binder of Stage 1
                                word for word) and singleTableLaw_of_moments
  General/Layer1.lean           the Layer 1 half of Stage 1, outside RMT/ (gate 8):
                                stack_law_general, prop_stacksvd_general_of_general,
                                lem_delocalization_of_general,
                                thm_svd_stack_general_of_general, thm_theta_est_of_general;
                                since Stage 3 also the same four with the edge supplied,
                                prop_stacksvd_general_of_moments,
                                lem_delocalization_of_moments,
                                thm_svd_stack_general_of_moments, thm_theta_est_of_moments
  RMT/Het/                      Gaussian discharge of HeteroLaw (route A, items H0 to H12):
                                MPhet, Split, Duality, R4het, R3het, Stein, R1het, R2het,
                                R5het, Simplicity, R6het, EdgeScalar, EdgeSharp
                                (heteroEdge_of_gaussian), Sup (heteroLaw_of_gaussian,
                                thm_stacksvd_weighted_gaussian, prop_dominance_gaussian)
  Secular.lean          lem:secular_equation, spectrum half (det_Rmat_sub, lamMax_Rmat_eq)
  ThetaEst.lean         thm:theta_est over SingleTableLaw and ThetaEstLaw (F31), item P
                        (noise_projection_tendsto), thetaEstLaw_of_gaussian; since
                        2026-09-08 also thetaEstLaw_of_general and thm_theta_est_general
                        (Stage 0: the two cross-table limits at NoiseLaw, the single-table
                        law of table i kept as a hypothesis)
  MLE.lean, MLEConverse.lean    app:wstacksvd_mle: the identity and its converse
  MLEMarginal/                  app:wstacksvd_mle, the Gaussian marginalization (L5, 2026-09-03):
                                Defs (reRowLaw, reTableLaw, reJointLaw, reDensity, sqrtMleCov),
                                RowLaw, TableLaw, LogLik, Main (thm_wstacksvd_mle_marginal)
  Remarks.lean, RemarksUniform.lean  the two remarks as explicit instances
  Sat.lean              satisfiability witnesses: Sat.Rank1.model, explicit Gaussian models
  Existence.lean        the existential forms: prop_binarystacksvd_inadmissable_exists,
                        remark_stack_outperform_svd_exists, remark_svd_outperform_stack_exists
  StrictFacades.lean    the strict clauses of prop:dominance and the binary comparison
  SVDStack/Rayleigh.lean        the bound uniform in w (paper finding E1)
  SVDStack/EntrywiseSigned.lean lem:entrywise_conv_eigenvec in the paper's signed form (L2)
  SVDStack/Inad.lean            prop:binarystacksvd_inadmissable, svdstack clauses (L3)
  SVDStack/Simple.lean          thm:simple_thm1, svdstack half, every M >= 1
  SVDStack/DelocFinite.lean     finite-d lem:delocalization tail bound 1/(eps^2(d-1)) (F1)
  LinAlg/KyFan.lean, Eigen.lean, Frame.lean     Ky Fan, padded spectra, top-r frames
  LinAlg/SpecIdx*.lean, SpecWindow.lean, SpecInvReindex.lean   index projector API
  RankR/                        Section 7 and Appendix E. Defs (r_i = 1 model), Unweighted
                                (its Layer 1 theorems, F28, 2026-09-08), General (general r_i
                                model, SpikedModelR, UnalignedModelR, TableLawR; since F28
                                also rk_le_d, norm_colVecG, the dedup survivors), GeneralMain,
                                GeneralFrob, GeneralGaussian, Weighted, WeightedMain (its
                                Layer 1 theorems, F28), WeightedUpper,
                                WeightedUpperG (thm:rank_r_svdstack aggregate), Subspace,
                                SubspaceMain (the Layer 1 theorems of Subspace and SubspaceG,
                                F28), SubspaceG, SubspaceGStack, SubspaceGaussian, Aligned,
                                AlignedComponent, AlignedMain (thm:rank_r_svdstack
                                component), AlignedOrth (svdstack columns orthogonal, F2),
                                StackGamma (HeteroLawR, since F28), StackMain
                                (thm:rank_r_stacksvd), Frobenius, Flatten, GramR, Example
  RankR/RMT/                    Gaussian TableLawR and SubspaceLawG at general r_i:
                                Stack, Split, Forms, Outliers, DelocR, DelocAffineR, EdgeR,
                                EdgeDetR, EdgeGlueR, EdgeGlueDetR, EdgeTauR, AlignG,
                                AlignOutG, AlignTauR, SimplicityR, SimplicityAffineR,
                                SymmetryR, DecompR, ShiftR, OutliersG, R6R, TableStack,
                                TableSimple, TableAlignSup, TableLawGaussian
  RankR/Het/                    Track E, the Gaussian discharge of HeteroLawR: Scalars,
                                Split, Duality, Simplicity, Edge (HeteroEdgeR), Forms
                                (ResolventLimitsHetR), Deloc, BulkDet, Outliers, Align,
                                Bulk, Sub (F8: the sub-model that drops the zero-weight
                                tables, UnalignedModelR.sub, key_of_gaussian), Sup
                                (heteroLawR_of_gaussian, thm_rank_r_stacksvd_gaussian)
  RankR/SingleWeight/           one weight per table (Appendix, 2026-09-05): Scalars (secMat,
                                IsSecularRoot, swTerm, swLimit, EigSep), ScalarsOne (the
                                r = 1 reduction), Defs (perfSW, vhatSW), Main
                                (SingleWeightLaw, prop_gen_rank_stacksvd_singleweight),
                                Example (swLimitEx, swLimitEx_lt), Existence
                                (witness_perfRG_tendsto,
                                prop_singleweight_suboptimality_of_law); Layer 1.
                                F18b (2026-09-07), the discharge of hlaw on the witness:
                                Regimes (EigSep on the witness where both roots are
                                detectable, swLimitEx = swLimit there, the tie value),
                                Tie (hlaw_tie: the weighted stack at equal weights is a
                                multiple of the unweighted one), Suboptimality (hlaw_both,
                                align_bulk_sw_one, hlaw_one, hlaw_witness,
                                prop_singleweight_suboptimality_gaussian, unconditional)
  RankR/SingleWeight/Het/       Track G, the Gaussian discharge of SingleWeightLaw at general
                                R_i (2026-09-05): Scalars (swFmat, swRho, swNu, the
                                reciprocal identity), Forms (ResolventLimitsSW, the column
                                and Gram limits), Frame, Count (deterministic: a frame of
                                approximate eigenvectors locates and counts the outliers),
                                Outliers (the count above tau and the projected norm),
                                Align (tendstoInProb_eigVal_sw, align_sw_of_gaussian), Sub
                                (F18c: the drop of the zero-weight tables, subRk, wSuppEmb,
                                simpleSpec_ae_stackGramW_of_exists, EigSep.exists_ne_zero),
                                Sup (singleWeightLaw_of_gaussian with
                                hrn : forall N, r <= sum over {i | w i /= 0} of n i N,
                                singleWeightLaw_of_gaussian_shift, the facades
                                prop_gen_rank_stacksvd_singleweight_gaussian and
                                _inner_gaussian; no condition on the weights, as in the
                                paper). F18b (2026-09-07), the one-detectable-root regime
                                of the witness: OneOutCount (at most one eigenvalue above
                                the edge), OneOutDet (deterministic: a positive secular
                                value caps a rank-one update, the second eigenvalue of a
                                two-column split, a form bound, Psi' -> +infinity),
                                OneOutAlign (the one outlier and its eigenvector),
                                OneOutBulk (the edge window at general R_i)
  Vendor/COLT83/                Stein, Sudakov-Fernique, Borell-TIS (backport of 2d84602)
```

Weyl, Davis-Kahan, Courant-Fischer and Gaussian concentration come from StatsMLlib, not from
us. The tree has no open proof since 2026-09-10 evening, when the last tracked count, the
Furedi-Komlos count, was proved (docs/SORRIES.md has no row); every declaration is proved
with the three standard axioms.

In the development repository, `CLAUDE.md` holds the working rules of the AI assistant and
`notes/SERVER_SETUP.md` the toolchain setup of the shared server and of the laptop; neither
is part of this copy.
