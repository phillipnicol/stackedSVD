# AUDIT_DOC.md: the entry point for the external audit

> Location (2026-09-03): this file lives in `docs/`. The bare names `TECHNICAL.md`, `THEOREMS.md`
> and `AUDIT_DOC.md` refer to the files of this folder; every other path is relative to the
> repository root. `paper_edits.md` and `NOTE_FOR_COAUTHORS.md` moved from this folder to
> `notes/` on 2026-09-11; this file names them `notes/paper_edits.md` and
> `notes/NOTE_FOR_COAUTHORS.md` throughout.
> The short `README.md` at the root is the plain-language entry point.

This file stands alone. Read it without the repository open and you learn what is claimed,
what is assumed, what is not proved, and how to check every claim about the Lean tree.

**Subject.** A Lean 4 formalization of *"Stacked SVD or SVD stacked? A random matrix theory
perspective on data integration"* (Baharav et al., AoS resubmission; arXiv:2507.22170 is the
older public version with the same labels) (`TECHNICAL.md`).

**How to read the sources.** Every factual line names the repository document it comes from, in
parentheses. Where two documents disagreed, both values were measured against the tree and the
stale one was corrected; no `CONFLICT` mark remains. The build and gate figures come from the
transcript in Appendix A, which is the verbatim output of one command.

---

## 1. The packet and how to read it

The packet has three files. This one is the entry point; the other two are annexes.

| File | What it is |
|---|---|
| `AUDIT_DOC.md` | this file: the claim, the assumptions, the gaps, and the check commands |
| `THEOREMS.md` | the per-result annex: paper statement (verbatim LaTeX), Lean signature (verbatim), hypotheses in words, modeling choices, scope note, audit trail, for every result; section 8 carries the verification record and section 8.1 the satisfiability witnesses |

**What is proved, in one sentence.** For `M` tables with independent Gaussian noise in the
proportional regime, the development gives machine-checked proofs of the Lean endpoints
listed in section 5, subject to the exact Lean model classes and definitions of `THEOREMS.md`
section 1, with the random matrix theory input proved rather than assumed, with no `axiom`
and no `sorry` (`TECHNICAL.md`, `THEOREMS.md` 2.7, `docs/SORRIES.md`). That is not the same as "the
paper's theorems": three paper results have no Lean statement (7.5), and two Lean statements
are not the paper's own (the general-`r_i` stackSVD result concludes on a projector sum and
not on `‖V̂ᵀ V‖_F²`, 7.8; the rank-`r` results assume a common component order in every table,
7.2 and 7.3). The rest of this document lists every such difference; the wording of this
paragraph follows the external statement audit of 2026-09-02
(`notes/archive/audit_external_statement_2026-09-02.md`, whose sign-off sentence it adopts, and
`notes/archive/audit_external_statement_response_2026-09-02.md`), with one change: that audit
listed a third difference, "the MLE result is an identity and not a derivation", which
closed on 2026-09-03 (7.15). Three of that audit's proposed additions are in the tree: the
entrywise eigenvector lemma in the paper's signed coordinate form (item L2, section 5, the
fourth difference of the earlier wording), the svdstack clauses of
`prop:binarystacksvd_inadmissable` on the model (item L3, 7.16), both since 2026-09-02, and
the Gaussian marginalization behind the MLE objective (item L5, 7.15, since 2026-09-03).

**Suggested reading order for an auditor.** Section 2 (the claim), section 3 (what is
assumed), section 7 (scope and gaps), then `THEOREMS.md` for the one or two results you want
to attack, then section 9 to run the gates yourself.

**What the packet does not contain.** The Lean sources. Every signature quoted in
`THEOREMS.md` is checked against the tree by a script, and section 9 says how to rerun that
check (`scripts/README.md`).

---

## 2. The claim

The claim is a set of named Lean theorems whose hypotheses are model hypotheses only:
independent Gaussian tables, unit-norm signal vectors, and the proportional regime
`n_i → ∞`, `d → ∞`, `n_i/d → c_i`. These are the `*_gaussian` names. They do not assume any
random matrix theory limit law: the limit laws that the paper cites as black boxes are proved
inside the development for Gaussian noise (`TECHNICAL.md`, "How the random matrix theory input
enters"; `THEOREMS.md` 2.7). Every conclusion is convergence in probability, never almost
surely (`THEOREMS.md` 0, 7.1).

Headline Gaussian endpoints, each with its paper label (`TECHNICAL.md` main table):

| Paper label | Gaussian theorem | File |
|---|---|---|
| `prop:single_table` | `singleTableLaw_of_gaussian` | `RMT/Full.lean` |
| `prop:stacksvd_general` | `prop_stacksvd_general_gaussian` | `StackSVD/Main.lean` |
| `thm:svd_stack_general` | `thm_svd_stack_general_gaussian` | `SVDStack/Main.lean` |
| `thm:stacksvd_weighted` | `thm_stacksvd_weighted_gaussian` | `RMT/Het/Sup.lean` |
| `thm:svdstack_weighted` | `thm_svdstack_weighted_gaussian` | `SVDStack/Weighted.lean` |
| `thm:theta_est` | `thm_theta_est_gaussian`; `thm_theta_est_general` (noise law `NoiseLaw`, 2026-09-08) | `ThetaEst.lean` |
| `prop:dominance` | `prop_dominance_gaussian` | `RMT/Het/Sup.lean` |
| `prop:stacksvd_subspace` | `prop_stacksvd_subspace_gaussian`, `_general_gaussian` | `RankR/SubspaceGaussian.lean` |
| `prop:general_rank_unweighted_svdstack` | `prop_general_rank_unweighted_svdstack_general_gaussian` | `RankR/GeneralGaussian.lean` |
| `thm:gen_rank_weight_svdstak` | `thm_gen_rank_weight_svdstak_general_r_full_gaussian` | `RankR/GeneralGaussian.lean` |
| `thm:rank_r_svdstack` | `thm_rank_r_svdstack_aggregate_gaussian`, `_component_gaussian` | `RankR/GeneralGaussian.lean`, `RankR/AlignedMain.lean` |
| `thm:rank_r_stacksvd` | `thm_rank_r_stacksvd_gaussian` | `RankR/Het/Sup.lean` |

Section 5 carries the mapping table, 26 rows, one per paper result that has a Lean statement,
with the Layer 1 name and the scope note of each. Every one of the 26 has a Gaussian name
since 2026-09-07 (`prop:singleweight_suboptimality`, `THEOREMS.md` 6.9, was the last: its
hypothesis `hlaw` is discharged by item F18b regime by regime, since the Gaussian discharge
of `SingleWeightLaw` alone does not supply it). The one paper result of 7.5, the removal
remark, has no row at all.

Two results are deterministic and carry no probability at all: `lem:secular_equation`
(`Secular.lean`, `StackSVDWeighted.lean`) and the MLE identity of `app:wstacksvd_mle`
(`MLE.lean`, `MLEConverse.lean`). A third, `lem:entrywise_conv_eigenvec`
(`SVDStack/Deterministic.lean`), has no Gaussian form because it needs none: its hypothesis is
entrywise convergence in probability of a random matrix sequence and its conclusion is a
`TendstoInProb` statement, so it is a statement about random matrices with no noise model
(`TECHNICAL.md` marks all three "deterministic" in the Gaussian column; for the third the word
means "no Gaussian discharge", corrected 2026-09-02).

---

## 3. What is assumed

### 3.1 Axioms: three, and a gate that enforces it

Every public declaration depends only on `propext`, `Classical.choice` and `Quot.sound`
(`TECHNICAL.md`). `CLAUDE.md` hard rule 2 forbids the `axiom` keyword anywhere in the
development, and `scripts/check_axioms.sh` is the gate that catches a new one
(`scripts/README.md`).

The gate is not a grep. `scripts/AxiomAudit.lean` imports the root module, walks every
declaration of the `StackedSVD` namespace, and calls `Lean.collectAxioms`, the function that
`#print axioms` itself uses (`scripts/README.md`). A declaration passes when its axiom set is
inside `{propext, Classical.choice, Quot.sound}`. `sorryAx` is the one conditional case: it
passes only when the declaration has a row in `docs/SORRIES.md`, or inherits one, and the driver
names the origin on a `SORRYDEPS` line (`scripts/README.md`).

### 3.2 No `sorry`

`docs/SORRIES.md` is empty since 2026-09-10 evening. Its last row was the Furedi-Komlos shape
count `cellCard_mul_le_pos`, now in `RMT/General/Edge/CountBound.lean` (`THEOREMS.md` 7.1;
empty before Stage 3, 2026-09-10), and `scripts/check_sorries.sh`
compares that ledger against the `sorry` sites in the tree. Run in this session:
`check_sorries: OK - 0 sorry sites in 0 declarations, 0 rows, all matched.`, exit 0. The
script blanks comments, doc comments and string literals before it matches, so the word
`sorry` in a comment is not a site (`scripts/README.md`).

Because `docs/SORRIES.md` is empty, a `sorryAx` anywhere in the tree now fails the axiom
gate (`THEOREMS.md` 7.1); no declaration carries one. Until 2026-09-10 evening the table
named exactly one declaration, the Furedi-Komlos count `cellCard_mul_le_pos`, and twelve
further declarations inherited `sorryAx` from it with no row of their own (the `SORRYDEPS`
mechanism, `scripts/README.md`), which the gate accepted while that row stood.

### 3.3 Vendored code cannot hide anything

A small set of Gaussian-comparison files is vendored from COLT83 (Apache 2.0): Stein's
identity, Gaussian interpolation, Sudakov-Fernique, Borell-TIS, 3716 lines in 20 files
(`TECHNICAL.md`, `notes/SERVER_TODO.md`; backport of commit `2d84602`, `CLAUDE.md`).
`check_sorries.sh` skips `Vendor/` and the axiom driver skips it unless
`AXIOM_AUDIT_VENDOR=1` (`scripts/README.md`), which invites the question: could a `sorry`
hide there?

It could not affect any result. `Lean.collectAxioms` is transitive, so a vendored lemma
carrying `sorryAx` or a new `axiom` would appear in the axiom set of every one of our
theorems that depends on it, and the gate would fail those theorems. Coverage of our results
therefore does not depend on auditing the vendored files directly. Direct greps also come back
empty: `grep -rn '\bsorry\b'` and the same grep for `axiom` over
`lean/StackedSVD/Vendor/` both return 0 (coordinator measurement, 2026-09-02).

### 3.4 Coverage: every module of ours is audited

The axiom driver only sees what `import StackedSVD` pulls in, so a module missing from the
root import file would be built and never audited (`scripts/README.md`).
`scripts/check_root_imports.py` closes that hole: it computes the transitive closure from
`lean/StackedSVD.lean` and fails when a module of ours is outside it.

On the build of record (commit `416ab9d` on `stage3`, the Stage 3 commit, 2026-09-10): 230 modules, 210 ours, 229 reachable, exit 0. The 197th to 210th modules of
ours are the fourteen modules of Stage 3 (`RMT/General/Edge/`, 2026-09-10: `Defs`, `Arith`,
`Trunc`, `Sparse`, `Trace`, `Compare`, `Gaussian`, `Count`, `Code` (the step code), `Dyck`
(the tree-walk bound), `CountBound` (the arithmetic assembly), `Excess`, `Markov`, `Sup`). The
one unreachable module is `Vendor/COLT83/Axioms.lean`, the vendor's own `#print axioms`
driver, which nothing imports on purpose (`scripts/README.md`). The 172nd to 197th modules of
ours were the 26 modules of the non-Gaussian Stage 1 (2026-09-09: `Prob/Chebyshev`,
`Prob/Tensorization`, 23 files under `RMT/General/`, `General/Layer1`); F35 retired
`RMT/General/Statements.lean` the same evening, which leaves 196. The 166th to 171st are the
six files of the import boundary F28 (`StackSVD/Main`, `StackSVD/Weighted`,
`SVDStack/DelocDir`, `RankR/Unweighted`, `RankR/WeightedMain`, `RankR/SubspaceMain`,
2026-09-08); the 163rd to 165th are the three `Prob/` files of Stage 0 (`NoiseLaw`,
`LinFormMoments`, `NoiseMoments`, the same day). Before that, at `8736b2a` (2026-09-08
morning): 182 modules, 162 ours, 181 reachable. The 162nd module of ours is `Main.lean` (2026-09-07
evening), the statement file, imported last by the root. The 155th to 161st modules of ours are the seven
F18b modules of `RankR/SingleWeight/` (`Regimes`, `Tie`, `Suboptimality`, `Het/OneOutCount`,
`Het/OneOutDet`, `Het/OneOutAlign`, `Het/OneOutBulk`, 2026-09-07); the 154th module of ours is
`RankR/SingleWeight/Het/Sub.lean` (F18c, 2026-09-05); the 147th to 153rd are the seven
other `RankR/SingleWeight/Het/` modules (Track G, the same day); the 146th is `RankR/Het/Sub.lean`
(F8, 2026-09-05); the 140th to 145th are the six `RankR/SingleWeight/` modules of the same
day; the 138th and 139th are `SVDStack/DelocFinite.lean` and `RankR/AlignedOrth.lean` (F1,
F2); the 137th is `Existence.lean` (section 6, 2026-09-03); the 127th was `Sat.lean`
(section 6); the 130th to 136th are the seven L5 modules of 2026-09-03 (7.15).

### 3.5 The noise model and the regime

One table is `X_N ω = θ u_N v_Nᵀ + d_N^{-1/2} Z_N ω` with `‖u_N‖ = ‖v_N‖ = 1`, `θ ≥ 0`, and
`Z` measurable (`StackedSVD/Defs.lean`, structure `SpikedModel`; `THEOREMS.md` 1.1).

- `gaussianMatrix n d` is the law of an `n × d` matrix with i.i.d. `N(0,1)` entries
  (`Defs.lean:129`).
- `SpikedModel.GaussianNoise` says `Z_N` has that law under `μ N` (`Defs.lean:168`).
- `MultiTableModel.JointGaussianNoise` says the joint law of the `M` noise matrices at level
  `N` is the product of `M` such laws (`Defs.lean:195`). This is the paper's
  `assum:general_noise` in the Gaussian case, plus the cross-table independence that the paper
  states in words (`THEOREMS.md` 1.1).
- `SpikedModel.Regime c` is a conjunction of three limits: `n → ∞`, `d → ∞`, and
  `n_N/d_N → c` (`Defs.lean:172`, `THEOREMS.md` 1.1). It is `eq:RMT_limit`.

There is one asymptotic index `N`, and one probability space `(Ω N, μ N)` per `N`
(`CLAUDE.md`, Conventions).

### 3.6 The model classes, in words

| Structure | File | What it fixes |
|---|---|---|
| `SpikedModel` | `Defs.lean` | one rank-one table: `θ ≥ 0`, unit `u_N`, unit `v_N`, measurable noise, `n_N > 0`, `d_N > 0` |
| `MultiTableModel` | `Defs.lean` | `M` tables that share the right singular vector `v` |
| `SpikedModelR` | `RankR/General.lean` | one table with `rk` spikes, orthonormal `U`, `V`; fields `hθnn : ∀ k, 0 ≤ θ k` (F8, 2026-09-05; `hθpos : ∀ k, 0 < θ k` before) and `hθanti : StrictAnti θ` |
| `UnalignedModel` | `RankR/Defs.lean` | `assum:unaligned` at `r_i = 1`: one unit vector `R_i` per table, tables built on `SpikedModel` |
| `UnalignedModelR` | `RankR/General.lean` | `assum:unaligned` at general `r_i`: `V` with `r` orthonormal columns, `R_i` with orthonormal columns, `V_i = V R_i` |

Two facts about this list matter for scope, and section 7 returns to both.
`SpikedModelR.hθnn` and `hθanti` are model fields, not theorem hypotheses; `hθanti` narrows
the class and `hθnn` (`0 ≤ θ`, F8) is wider than the paper's positive entries (`THEOREMS.md`
1.1). `UnalignedModel` is not `UnalignedModelR` at `rk = 1`, and the
tree carries no bridge between them, so a general-`r_i` theorem does not visibly imply its
`r_i = 1` twin; both twins are proved separately (`THEOREMS.md` 1.1).

### 3.7 The random matrix theory input is proved, not assumed

**This is the sentence an auditor most needs.** The paper cites `liu2023asymptotic` Theorems 1
and 2 and `10.3150/19-BEJ1129` Theorem 2.3 as black boxes. No prover reaches those today, so
they enter this development as `structure`s, one field per fact a proof reads, and each
structure is then **proved** for Gaussian noise by a named theorem. The `*_gaussian` endpoints
of section 2 therefore carry no limit-law hypothesis (`THEOREMS.md` 2, 2.7; `TECHNICAL.md`).

The facades, with their files (`THEOREMS.md` 2.7):

| Structure | Defined in | Gaussian discharge | Discharge file |
|---|---|---|---|
| `SpikedModel.SingleTableLaw` | `RMT.lean` | `singleTableLaw_of_gaussian` | `RMT/Full.lean` |
| `MultiTableModel.HeteroLaw` | `StackSVDWeighted.lean` | `heteroLaw_of_gaussian` | `RMT/Het/Sup.lean` |
| `MultiTableModel.HeteroEdge` | `RMT/Het/R4het.lean` | `heteroEdge_of_gaussian` | `RMT/Het/EdgeSharp.lean` |
| `SpikedModelR.TableLawR` | `RankR/General.lean` | `tableLawR_of_gaussian_rk` | `RankR/RMT/TableLawGaussian.lean` |
| `UnalignedModel.SubspaceLaw` | `RankR/Subspace.lean` | `subspaceLaw_of_gaussian` | `RankR/SubspaceGaussian.lean` |
| `UnalignedModelR.SubspaceLawG` | `RankR/SubspaceG.lean` | `subspaceLawG_of_gaussian` | `RankR/SubspaceGaussian.lean` |
| `UnalignedModelR.HeteroLawR` | `RankR/StackMain.lean` | `heteroLawR_of_gaussian` | `RankR/Het/Sup.lean` |
| `UnalignedModelR.HeteroEdgeR` | `RankR/Het/Edge.lean` | `heteroEdgeR_of_gaussian` | `RankR/Het/Edge.lean` |
| `UnalignedModelR.ResolventLimitsHetR` | `RankR/Het/Forms.lean` | `resolventLimitsHetR_of_gaussian` | `RankR/Het/Forms.lean` |
| `UnalignedModelR.SingleWeightLaw` | `RankR/SingleWeight/Main.lean` | `singleWeightLaw_of_gaussian` (under the paper's separation assumption `EigSep`, a condition on the parameters, and `hrn`) | `RankR/SingleWeight/Het/Sup.lean` |

Every limit-law hypothesis structure in the tree has a Gaussian discharge (`EigSep` is a
condition on the parameters, not a law, and stays a binder), and no result of the paper
that this development covers rests on an undischarged black box (`THEOREMS.md` 2.7).

Two of these were expected to be out of reach. The exact heteroscedastic upper edge, which two
expert reviews put beyond Bai-Silverstein, is proved by a sharp Sudakov-Fernique bound
(`heteroEdge_of_gaussian`, `THEOREMS.md` 2.2). `heteroLaw_of_gaussian_margin` is the weaker
interim theorem and is kept only as a record; it needs `MPhet.bSF c w < MPhet.rhoHet θ c w`
(`THEOREMS.md` 2.2).

### 3.8 Convergence in probability is the paper's own mode, not a downgrade

Every Lean conclusion is `TendstoInProb`, convergence in probability
(`THEOREMS.md` 0, 1.2). The paper claims the same mode, so the formal claim is the paper's
claim. In the read-only snapshot `main_paper.tex`: line 234 defines `\pto` as "convergence in
probability of a sequence of random variables"; line 159 says "We establish convergence in
probability to the limiting forms"; line 175 repeats it for the rank-one section. The words
"almost surely" occur twice in the paper, at lines 1129 and 1142, and both are inside one
proof, not in a theorem statement (coordinator verification, 2026-09-02; the three lines and
the two occurrences were re-read in the snapshot this session).

Nothing in the development is stated as an almost-sure limit in `N`. Two facts inside the
hypothesis structures are almost sure at each fixed `N` (`SingleTableLaw.topSimple`,
`TableLawR.simple`); those are per-`N` statements and strengthenings, not limits
(`THEOREMS.md` 7.1).

---

## 4. The two-layer design, and why it is not circular

The design has two layers (`notes/archive/PLAN_2026-08-29.md`, `THEOREMS.md` 0):

1. **Layer 1** derives the paper's conclusion from a hypothesis structure. It is an honest
   implication and holds for any noise law that satisfies the structure.
2. **Layer 2** proves the structure for Gaussian noise. Its name ends in `_of_gaussian`.

A `*_gaussian` endpoint is Layer 1 applied to Layer 2. Nothing is left assumed.

**Dependency direction, which is where a reader checks for circularity.** For the homogeneous
chain the file layout shows it. `RMT.lean`, which defines `SingleTableLaw`, imports only
`StackedSVD.Defs`; `RMT/Full.lean`, which proves `singleTableLaw_of_gaussian`, imports
`RMT.TailShift` and `RMT.R6` (the Marchenko-Pastur chain) and no Layer 1 file;
`StackSVD/Main.lean` (through `StackSVD.lean`) and `SVDStack/Main.lean`, which state the
paper results, import `RMT.Full` (import headers read 2026-09-03, `StackSVD/Main.lean`
since F28, 2026-09-08). The Marchenko-Pastur chain that proves the single-table law is `RMT/{MP,
MP7, R0, R4, R4C, Simplicity, Symmetry, R3, R3minus, T, R1, ResolvDeriv, SteinStep, R2, R5,
R6, Sup, TailShift}`, in that order (`THEOREMS.md` 2.1).

For the heteroscedastic chain the file layout shows the split less directly. Before F28
(2026-09-08) `StackSVDWeighted.lean` held the weighted definitions, the structure
`HeteroLaw` and the Layer 1 consumer `thm_stacksvd_weighted_general` in one file; F28 moved
the consumer, unchanged, into `StackSVD/Weighted.lean`, so `StackSVDWeighted.lean` now holds
only the definitions and `HeteroLaw`. Three files of the discharge still import
`StackSVDWeighted.lean` for those definitions: `RMT/Het/MPhet.lean`, `RMT/Het/Split.lean`
and `RMT/Het/Simplicity.lean` (external rank-one audit of 2026-09-03, finding 10). A fourth
file, `RMT/Het/Sup.lean`, imports `StackSVD/Weighted.lean` as well, since it proves the
final Gaussian corollary, `thm_stacksvd_weighted_gaussian`, which needs the Layer 1 theorem
by design; that file is meant to mix a discharge and a consumer, and it is one of three
such files gate 8 (`scripts/check_core_imports.py`) keeps outside its protected set, EXT
(section 9). Lean forbids a declaration cycle, so `heteroLaw_of_gaussian`
cannot use a theorem that uses it; what the import does not rule out is a discharge that
reuses an earlier Layer 1 theorem about the same structure. The declaration-level check
`scripts/check_layering.sh` (driver `scripts/LayerAudit.lean`) rules that out: it computes
the transitive constant closure of each of the three discharges (`singleTableLaw_of_gaussian`,
`heteroLaw_of_gaussian`, `heteroEdge_of_gaussian`) and fails if the closure contains any
declaration whose type takes a `SingleTableLaw` or a `HeteroLaw` as a hypothesis. It prints
the closure sizes and the `USES` records (which Layer 1 constants each discharge does touch).
It also flags 13 Layer 1 theorems by name, and self-tests the type criterion on five known
consumers before the audit. Result on the build of record: 3 endpoints audited, closures of 1288 (`singleTableLaw_of_gaussian`), 1562 (`heteroLaw_of_gaussian`) and 323 (`heteroEdge_of_gaussian`) project constants, 97 `USES` lines over 75 distinct Layer 1 constants (the structure `HeteroLaw`, its constructor, and definitions and structural lemmas of `StackSVDWeighted.lean`; no theorem that takes a structure as a hypothesis), 0 `FORBIDDEN`, 0 `CONSUMER`, exit 0; the counts are those of the build of record of 2026-09-05 (19:24 to 19:26 EDT; 98 `USES` lines since the F32 dedup of 2026-09-08, see section 9), the 30 s timing is from the run of 2026-09-03 09:10 EDT.
The heteroscedastic chain is `RMT/Het/` (items H0 to H12) and its rank-`r` twin is
`RankR/Het/` (Track E, stages E0 to E10) (`THEOREMS.md` 0 glossary, 2.5). The file split of
`StackSVDWeighted.lean` (F28, done 2026-09-08) saved no lines and the declaration-level
check stays the stronger statement regardless of file layout.

The mechanical check that closes the `sorry` question is the axiom gate, not the file layout:
`collectAxioms` is transitive, so a Gaussian endpoint that still depended on an unproved
statement would carry `sorryAx` and fail (`scripts/README.md`, `THEOREMS.md` 8). The axiom
gate does not check the layering; `check_layering.sh` does.

**One caution about the phrase "Layer 1".** A theorem is noise-free only when its signature
carries no noise hypothesis. The whole SVDstack family (`thm_svd_stack_general`,
`thm_svdstack_weighted` and their variants) takes both a `SingleTableLaw` and
`hI : m.IndepNoise`, because the cross-table step `lem_delocalization` needs a Fubini step
over the product law of the tables. `IndepNoise` asks only for independence across tables
under some product law with arbitrary marginals (`SVDStack/Gram.lean`, since 2026-09-02, L4;
before that the binder was `hG : m.JointGaussianNoise`); at rank `r` the Track B statements
take the same predicate as `hG : m.IndepNoise`. That is weaker than Gaussian, and it is not
the same hypothesis as `JointGaussianNoise` (`notes/archive/AUDIT_PACK_V4.md` 2.4, C17).
`thm_theta_est` reads its two cross-table limits through `MultiTableModel.ThetaEstLaw`
(`THEOREMS.md` 2.6) since F31, 2026-09-08; until then it was the one rank-one Layer 1 theorem
that kept `hG : m.JointGaussianNoise` (section 7.1).

---

## 5. Paper result to Lean theorem

The table below maps every paper result that has a Lean statement, 26 rows; the one result of
7.5 has no row (`TECHNICAL.md`). Column meanings: **Gaussian theorem** takes model hypotheses
only, or names the Layer 1 discharge that is still open; **Layer 1 theorem** takes the random
matrix theory input as an explicit hypothesis structure; **File** names the file of the main
declaration of the row, so a scalar half in `Scalars.lean` or a Gaussian discharge in `RMT/`
or `RankR/RMT/` may live elsewhere (`THEOREMS.md` gives the file of every declaration);
**Scope note** is the one-line summary of the difference from the paper.

Reading the names: a name that starts with `_` continues the name before it. In the row for
`thm:svd_stack_general` the Layer 1 column reads `thm_svd_stack_general`, `_inner`, `_zero`,
that is three declarations `thm_svd_stack_general`, `thm_svd_stack_general_inner` and
`thm_svd_stack_general_zero` (`TECHNICAL.md`).

| Paper label | Gaussian theorem | Layer 1 theorem | File | Scope note |
|---|---|---|---|---|
| `prop:single_table` | `singleTableLaw_of_gaussian` | `SpikedModel.SingleTableLaw` (structure) | `RMT/Full.lean`, `RMT.lean` | both regimes; `topSimple` added |
| `lem:delocalization` | via `lem_delocalization` + `singleTableLaw_of_gaussian` (`RMT/Full.lean`) | `lem_delocalization` | `SVDStack/Gram.lean` | signed limit `⟪v̂_i, v̂_j⟫ → β_iβ_j` |
| `lem:entrywise_conv_eigenvec` | none needed (the hypothesis is entrywise convergence in probability) | `lem_entrywise_conv_eigenvec`, `lem_entrywise_conv_eigenvec_signed` | `SVDStack/Deterministic.lean`, `SVDStack/EntrywiseSigned.lean` | projector form, sign-free, and the paper's signed coordinate form (added 2026-09-02) |
| `prop:stacksvd_general` | `prop_stacksvd_general_gaussian`, `_inner_gaussian` | `prop_stacksvd_general`, `_inner` | `StackSVD/Main.lean` | reduces to the stack at `(‖θ‖₂, ‖c‖₁)` |
| `thm:svd_stack_general` | `thm_svd_stack_general_gaussian`, `thm_svd_stack_general_inner_gaussian`, `thm_svd_stack_general_zero_gaussian` | `thm_svd_stack_general`, `_inner`, `_zero` | `SVDStack/Main.lean` | `β_2 > 0` unsorted |
| `thm:simple_thm1` (cor. 1) | `thm_simple_thm1_stacksvd_gaussian`, `thm_simple_thm1_svdstack_gaussian_full` | `Scalars.simple_thm1_stacksvd`, `_svdstack` | `StackSVD/Main.lean`, `SVDStack/Simple.lean` | svdstack half for every `M ≥ 1` |
| `cor.2` (binary) | `stackPerfW_binary_tendsto_gaussian`, `exists_binary_tendsto_max_gaussian` | `Scalars.binaryStackSVDLimit_eq_stackSVDLimit` | `StackSVD/Weighted.lean` | any nonempty subset, strict threshold |
| `thm:stacksvd_weighted` | `thm_stacksvd_weighted_gaussian`, `_gaussian_opt`, `_gaussian_inner`, `_gaussian_smul` | `thm_stacksvd_weighted`, `_inner`, `_general` | `RMT/Het/Sup.lean`, `StackSVD/Weighted.lean` | both regimes, exact edge, `∃ i, θ_i ≠ 0` |
| `thm:svdstack_weighted` | `thm_svdstack_weighted_gaussian`, `_gaussian_opt`, `thm_svdstack_weighted_gaussian_opt_full`, `thm_svdstack_weighted_paper_gaussian` | `thm_svdstack_weighted`, `_paper`, `_general`, `_inner`, `_zero` | `SVDStack/Weighted.lean`, `SVDStack/Rayleigh.lean` | `_opt_full` is uniform in `w` (finding E1) |
| `lem:secular_equation` | deterministic | `secular`, `IsGammaTop`, `det_Rmat_sub`, `lamMax_Rmat_eq`, `topSimple_Rmat` | `StackSVDWeighted.lean`, `Secular.lean` | finite matrix algebra |
| `thm:stacksvd_binary_optimal_svd_stack` | `thm_stacksvd_binary_optimal_svd_stack_gaussian`, `thm_stacksvd_binary_optimal_svd_stack_gaussian_strict` | `Scalars.svdstackOpt_le_binary`, `Scalars.svdstackOpt_lt_binary` | `StackSVD/Weighted.lean`, `StrictFacades.lean` | `c_i ≤ 1` on kept tables only |
| `prop:dominance` | `prop_dominance_gaussian`, `prop_dominance_gaussian_strict`, `_strict_svdstack`, `_strict_unweighted` | `Scalars.svdstackOpt_le_stackSVDLimitW`, `stackSVDLimit_le_stackSVDLimitW` | `RMT/Het/Sup.lean`, `StrictFacades.lean` | both dominance clauses |
| `prop:binarystacksvd_inadmissable` | `prop_binarystacksvd_inadmissable_exists` (the paper's existential form, model built), `prop_binarystacksvd_inadmissable_gaussian` (the five limits on any model of the instance), `stackPerfW_binary_inad_gaussian`, `stackPerfW_opt_inad_gaussian`, `svdstackPerf_inad_gaussian`, `svdstackPerfW_opt_inad_gaussian` | `Scalars.binarystacksvd_inadmissable`, `_ceil` | `Existence.lean`, `RMT/Het/Sup.lean`, `SVDStack/Inad.lean`, `Scalars.lean` | explicit instance, uniform over subsets; the two svdstack clauses on the model added 2026-09-02; the existential form with the model `Sat.Rank1.inadModel M` (`k_i = 2i + 1`, `l = 1`) added 2026-09-03 |
| `thm:theta_est` | `thm_theta_est_gaussian`; `thm_theta_est_general` (any fixed noise law with mean 0, variance 1, finite fourth moment; the single-table law of table `i` stays a hypothesis, 2026-09-08) | `thm_theta_est` | `ThetaEst.lean` | item P proved, not assumed; `ThetaEstLaw` also discharged at `NoiseLaw` (`thetaEstLaw_of_general`) |
| `app:wstacksvd_mle` | deterministic | `mleLogLik_eq`, `mleLogLik_max_iff_mem_topSpace`, `thm_wstacksvd_mle_marginal` | `MLE.lean`, `MLEConverse.lean`, `MLEMarginal/Main.lean` | finite-`N` identity, its converse, and the Gaussian marginalization from the random-effects model (L5, 2026-09-03) |
| `remark:stack_outperform_svd` | `remark_stack_outperform_svd_exists` (model built, every `M ≥ 2`), `remark_stack_outperform_svd_stack`, `_svdstack`, `_svdstack_uniform` | instances of the closed forms | `Existence.lean`, `Remarks.lean`, `RemarksUniform.lean` | explicit rationals; the existential form with the model `Sat.Rank1.oneModel M` added 2026-09-03 |
| `remark:svd_outperform_stack` | `remark_svd_outperform_stack_exists` (model built, `θ = (√5, 4)`, `c = (1, 38.4)`), `remark_svd_outperform_stack_two`, `_pair_*`, `_three_*` | instances of the closed forms | `Existence.lean`, `Remarks.lean`, `RemarksUniform.lean` | the `M = 3` example is on the threshold (finding E6); the existential form with the model `Sat.Rank1.twoModel` (`n = (5(N+1), 192(N+1))`, `d = 5(N+1)`) added 2026-09-03 |
| `lem:general_rank_delocalization` | via `tableLawR_of_gaussian_rk` (`RankR/RMT/TableLawGaussian.lean`) | `lem_general_rank_delocalization_general`, `lem_general_rank_delocalization` | `RankR/GeneralMain.lean`, `RankR/Defs.lean` | off-diagonal half; `IndepNoise` suffices |
| `prop:general_rank_unweighted_svdstack` | `prop_general_rank_unweighted_svdstack_general_gaussian`, `prop_general_rank_unweighted_svdstack_gaussian`, `prop_general_rank_unweighted_svdstack_general_frobenius_eig_gaussian` | `prop_general_rank_unweighted_svdstack_general`, `_frobenius`, `_frobenius_eig` | `RankR/GeneralGaussian.lean`, `RankR/GeneralFrob.lean`; the `r_i = 1` twin in `RankR/Defs.lean` | general `r_i`; trace form and Frobenius form |
| `prop:stacksvd_subspace` | `prop_stacksvd_subspace_gaussian`, `prop_stacksvd_subspace_general_gaussian`, `prop_stacksvd_subspace_frobenius_eig_gaussian` | `prop_stacksvd_subspace`, `_general`, `_frobenius`, `_frobenius_eig` | `RankR/SubspaceGaussian.lean`, `RankR/Subspace.lean` | no eigengap of `C`; ordered spikes per table |
| `thm:gen_rank_weight_svdstak` | `thm_gen_rank_weight_svdstak_general_r_paper_gaussian`, `thm_gen_rank_weight_svdstak_general_r_full_gaussian`, `thm_gen_rank_weight_svdstak_gaussian` | `thm_gen_rank_weight_svdstak`, `_opt`, `_max`, `_paper`, `_general_r`, `_full` | `RankR/GeneralGaussian.lean`, `RankR/Weighted.lean` | `rank B_R = r` in place of `β_ij > 0`; uniform-in-`W` bound (finding E2) |
| `eq:psi_equation`, `eq:perf_rankr_ex_*` | `example_perfR_tendsto_gaussian`, `example_perfStackR_tendsto_gaussian` | `example_perfR_tendsto` (`RankR/Example.lean`), `example_perfStackR_tendsto` (`RankR/Subspace.lean`) | `RankR/Example.lean`, `RankR/SubspaceGaussian.lean` | the worked example and its kinks |
| `thm:rank_r_svdstack` | `thm_rank_r_svdstack_aggregate_gaussian`, `thm_rank_r_svdstack_component_gaussian`, `_of_model_gaussian` | `thm_rank_r_svdstack_aggregate`, `_component`, `_component_eig` | `RankR/GeneralGaussian.lean`, `RankR/WeightedUpperG.lean`, `RankR/AlignedMain.lean` | both clauses; canonical frame; not a corollary (finding E4) |
| `thm:rank_r_stacksvd` | `thm_rank_r_stacksvd_gaussian`, `thm_rank_r_stacksvd_proj_gaussian`, `_inner_gaussian`, `_frobenius_gaussian` | `thm_rank_r_stacksvd`, `_proj`, `_inner`, `_frobenius` | `RankR/Het/Sup.lean`, `RankR/StackMain.lean` | ordered spikes, `R_i = 1`, so `ℓ_j = j` (finding E7) |
| `prop:gen_rank_stacksvd_singleweight` | `prop_gen_rank_stacksvd_singleweight_gaussian`, `_inner_gaussian`; `singleWeightLaw_of_gaussian` (with `hrn : ∀ N, r ≤ ∑ i with w i ≠ 0, n i N`), `singleWeightLaw_of_gaussian_shift` (without it) | `prop_gen_rank_stacksvd_singleweight`, `_inner`; `EigSep`, `SingleWeightLaw` (structures); `exists_unit_eigvec_secMat` (deterministic first half) | `RankR/SingleWeight/Het/Sup.lean` (and the seven files before it in `Het/`), `RankR/SingleWeight/Scalars.lean`, `Defs.lean`, `Main.lean` | one weight `w_i` per table, general `R_i`; matrix secular equation on `secMat`; `EigSep.sorted` excludes ties (finding E11); no condition on the weights, as in the paper (F18c dropped D37's `w_i ≠ 0` on 2026-09-05; `Het/Sub.lean` drops the zero-weight tables first) |
| `prop:singleweight_suboptimality` | `prop_singleweight_suboptimality_gaussian` (F18b, 2026-09-07: `hlaw_witness` proves the single-weight limit on the witness for every positive weight pair, the tie by the unweighted rank-`r` law, the two-root pairs by Track G, the one-root pairs by a one-component outlier law and an edge window) | `prop_singleweight_suboptimality_of_law`, with `hlaw` the one hypothesis | `RankR/SingleWeight/Example.lean`, `Existence.lean`, `Suboptimality.lean` | existence instance `M = 2`, `r = 2`, `θ_0 = 8/5`; every `w_i > 0` stays strictly below unweighted svdstack's `2 β_0²`; exact limit on the one-root pairs (E12) |

Naming conventions used above (`THEOREMS.md` 0 glossary): a **facade** restates a proved
result in the paper's own shape (inner product in place of projector, one weighting in place
of another) and adds no mathematics; facade names end in `_inner`, `_paper`, `_zero`, `_opt`,
`_full`, `_frobenius`, `_eig`. Overlaps are stated in projector form, which is sign-free; the
`_inner` twin is the paper's own display and consumes a simplicity fact
(`THEOREMS.md` 1.2, `TECHNICAL.md`).

---

## 6. The hypothesis sets are satisfiable

A theorem whose hypotheses no model meets is true and empty. `lean/StackedSVD/Sat.lean`
answers that objection inside the tree: it builds concrete Gaussian models and applies a
Gaussian endpoint to each, with no hypothesis left over (`Sat.lean` module docstring;
`THEOREMS.md` 8.1; `TECHNICAL.md`). None of the three headline endpoints is vacuous.

| Witness theorem | Endpoint applied | Paper label |
|---|---|---|
| `StackedSVD.Sat.Rank1.sat_stacksvd_weighted` | `thm_stacksvd_weighted_gaussian` | `thm:stacksvd_weighted` |
| `StackedSVD.Sat.Rank1.sat_svdstack_weighted` | `thm_svdstack_weighted_gaussian` | `thm:svdstack_weighted` |
| `StackedSVD.Sat.RankR.sat_rank_r_stacksvd` | `thm_rank_r_stacksvd_gaussian` | `thm:rank_r_stacksvd` |

Since 2026-09-03 the same construction supplies the models of the paper's three existential
statements (`Existence.lean`, external rank-one audit of 2026-09-03, findings 1 and 2):
`Sat.Rank1.inadModel M` (`k_i = 2i + 1`, `l = 1`) for `prop_binarystacksvd_inadmissable_exists`,
`Sat.Rank1.oneModel M` (`k_i = l = 1`) for `remark_stack_outperform_svd_exists`, and
`Sat.Rank1.twoModel` (`k = (5, 192)`, `l = 5`, so `c = (1, 38.4)`) for
`remark_svd_outperform_stack_exists`. Each of these theorems concludes `∃ ... (m :
MultiTableModel μ M n d), ...` with no hypothesis on a model (section 5, `THEOREMS.md` 5.3
and 5.6).

**The rank-one construction is parameterized**, not a single hand-tuned instance: `model k l
hk hl θ hθ` builds `M` tables with `n_i N = k_i (N+1)` and `d N = l (N+1)`, so the aspect ratio
is exactly `c_i = k_i / l` at every `N`, with noise `Measure.pi (fun i => gaussianMatrix (n i
N) (d N))` (`Sat.lean`). The instance used is `M = 2`, `c = (1, 1)`, `θ = (2, 2)`, so
`θ_i⁴/c_i = 16 > 1` and both tables are above the detection threshold (`Sat.lean`).

**The rank-`r` witness** is `M = 2`, `r = 2`, `n_i N = d N = N + 3` (so `c_i = 1`),
`θ = [[3, 3/2], [5/2, 1]]`, `R_i = 1` (`Sat.lean`). Both components are supercritical
(`∑_i θ_i0⁴/c_i = 120.06`, `∑_i θ_i1⁴/c_i = 6.06`), so both limits are strictly positive:
`γ_0 = 0.928594873932` and `γ_1 = 0.614906970741`, matched to `1.0e-13` against
`compute_x_star` of the R package `theory_pred.R`
(`Sat.lean`, `notes/archive/audit_independent_TrackE_wave2b_2026-09-02.md` sections 5.3 and 5.4).

Three facts about this file:

1. It compiles: commit `1b56d91`, `lean-local.sh`, exit 0, no errors (coordinator
   measurement, 2026-09-02).
2. It is inside the import closure of `lean/StackedSVD.lean` (line 129), and so is
   `Existence.lean` (line 130), so `check_axioms.sh` audits the witnesses and the three
   existential theorems on every run and they cannot rot (`Sat.lean` docstring;
   `check_root_imports.py`).
3. The models come from two independent audits that built them outside the tree: the vacuity
   probe of `notes/archive/audit_independent_D33_2026-09-02.md` and section 5.3 of
   `notes/archive/audit_independent_TrackE_wave2b_2026-09-02.md` (`Sat.lean` docstring).

Scope of a witness: it settles satisfiability of one hypothesis set. It says nothing about
other models and it changes no statement (`Sat.lean` docstring).

**One hypothesis structure has no witness in `Sat.lean`.** `SingleWeightLaw` (`THEOREMS.md`
6.8) is not in `Sat.lean`. The natural candidate is a rank-one model with one table at
`θ = 0`, and no `MultiTableModel`-to-`UnalignedModelR` bridge exists yet to build it (since
F8, 2026-09-05, `SpikedModelR.hθnn` allows `θ = 0` and `SpikedModel.toRankR` maps one
rank-one table; the multi-table bridge is still to write). Since Track G (2026-09-05) the
Gaussian discharge `singleWeightLaw_of_gaussian` builds the structure on every Gaussian
`UnalignedModelR` with `r ≤ ∑ i with w i ≠ 0, n i N` at every `N` that meets `EigSep`, so the hypothesis set is
consistent wherever such a model exists; the two-table witness `mdl` of `THEOREMS.md` 6.9
at `w = (1, 9/10)` meets `EigSep` numerically (`THEOREMS.md` 8.1), not yet in Lean.

---

## 7. Scope and known gaps

This section is complete to the best of the repository's own record. Each item leads with the
limitation.

**7.1 Gaussian noise only. This is the largest single gap.** Every `*_gaussian` theorem
assumes Gaussian entries. The paper's `assum:general_noise` asks only for i.i.d. entries with
mean 0, `E(√d z)² = 1` and a bounded fourth moment (`notes/README.md`). Universality
from Gaussian noise to that model is not formalized (`notes/archive/AUDIT_PACK_V4.md` S1). The
mitigation is structural, not proved, and it does not cover every family. The stackSVD Layer
1 theorems (rank one and rank `r`) and the rank-`r` SVDstack theorems read the noise only
through their limit-law structure and, at rank `r`, through `UnalignedModelR.IndepNoise`
(independent tables with arbitrary laws, `RankR/GramR.lean`); a four-moment discharge of the
same structures would give those results with no change to any statement. The rank-one
SVDstack family (`thm_svd_stack_general`, `thm_svdstack_weighted`, the uniform bound of
`SVDStack/Rayleigh.lean`, and `lem_delocalization` under them) reads the noise through
`MultiTableModel.IndepNoise` (`SVDStack/Gram.lean`) since 2026-09-02 (L4 of the external
statement audit, `notes/archive/L4_indepnoise.md`): the Fubini step
`measure_deloc_le_of_pi_indep` is stated for an arbitrary product law, so the same holds
there. `thm_theta_est` reads the noise through `MultiTableModel.ThetaEstLaw` (`THEOREMS.md`
2.6) since F31, 2026-09-08: item P with a random direction (`‖E_j v̂_i‖² → c_j`) and the
cross term (`u_jᵀ E_j v̂_i → 0`), proved for Gaussian noise as Chebyshev bounds with second
moments under `gaussianMatrix` (`noise_projection_topDir_tendsto`, `cross_term_tendsto`);
with arbitrary marginals both limits fail, so a non-Gaussian discharge of `ThetaEstLaw`
needs a moment hypothesis on the marginals, not only independence (`THEOREMS.md` 7.1,
`notes/FLAGGED.md` D35). Until F31 the Layer 1 form kept `hG : m.JointGaussianNoise` itself,
by user decision Q3 of 2026-08-29 (`notes/FLAGGED.md` Q3).

**7.2 The rank-`r` model class orders the spikes inside a table, and the paper does not.**
`SpikedModelR.hθnn : ∀ k, 0 ≤ θ k` (F8) and `hθanti : StrictAnti θ` are model fields. The paper
makes no ordering assumption (`main_paper.tex:767`) and assumes distinct entries for svdstack
only (`:768`); it says `prop:stacksvd_subspace` applies with repeated eigenvalues (`:842`). So
`prop_stacksvd_subspace_general_gaussian` is narrower than the paper at a table with a
repeated `θ_ij`. No stacksvd proof reads either field, so an order-free model would carry the
same proof (`THEOREMS.md` 6.3, Track A audit finding F1, decision D-A1 in `notes/FLAGGED.md`
item 10).

**7.3 Inside the ordered class the index `ℓ_j` collapses to `j`, and the paper's own worked
example is outside the class.** The spectral index of `thm:rank_r_stacksvd` is
`Scalars.ellR m.thetaAligned c j`, the paper's `ℓ_j`. It equals `j` only through
`ellR_thetaAligned` (`RankR/StackGamma.lean:316`), which takes **two side conditions**,
`hc : ∀ i, 0 < c i` and `hM : 0 < M`. The Layer 1 theorem `thm_rank_r_stacksvd` assumes
neither, so read its index as the paper's `ℓ_j`; the Gaussian form
`thm_rank_r_stacksvd_gaussian` has both, so there `ℓ_j = j` (`THEOREMS.md` 6.7). The paper's
worked example of `alg:rank_r_stacksvd` (`main_paper.tex:2359`), `θ_1 = (2, 1)` and
`θ_2 = (1, 10)`, is therefore **not** an `UnalignedModelR` with `R_i = 1` (`THEOREMS.md` 6.7).
Whether to keep the class or add a permutation `R_i` route is open question Q16 item 1
(`notes/FLAGGED.md`). The same common order narrows `thm:rank_r_svdstack`:
`thm_rank_r_svdstack_aggregate_gaussian` and `_component_gaussian`
(`RankR/GeneralGaussian.lean:201, 219`) take the same `UnalignedModelR μ M n d r (alignedRk M r)`
with `hR : ∀ i, m.R i = 1`, so every table ranks the `r` global components in the same strict
order there too; the hypothesis `hS : StrictAnti (Sagg β)` is a relabeling, the per-table
`hθanti` is the restriction (added 2026-09-02, external statement audit finding 7).

**7.4 At unordered `θ` the paper's first display is false, which is finding E7.** With the
paper's own `ℓ_j` (the rank of `θ̃_jj`), a Monte Carlo instance at `M = 2`, `r = 2`,
`c = (1, 1)` puts `0.571 ± 0.008` at index 0 and `0.011 ± 0.003` at the paper's index 1 at
`d = 1200`. Inside the ordered class the display holds (`THEOREMS.md` 6.7,
`notes/paper_edits.md` E7, `notes/archive/audit_rankr_plan_E_2026-09-02.md` row 4). The evidence is a
closed-form computation of the outlier locations plus Monte Carlo for the index assignment
and the limits; no Lean theorem covers unordered `θ`, so E7 is a numerically supported
counterexample and not a certified one (external statement audit finding 14).

**7.5 One result of the paper has no Lean statement; two more are formalized at Layer 1
only** (`THEOREMS.md` 7.1):

| Paper result | Label | Why not |
|---|---|---|
| The removal remark | `main_paper.tex:901` | not formalized, and wrong in general (Q8) |

Two more results, `prop:singleweight_suboptimality` (`main_paper.tex:915`) and
`prop:gen_rank_stacksvd_singleweight` (`:2112`), have Lean statements since 2026-09-05
(`THEOREMS.md` 6.8, 6.9, section 5 above). The second is proved at Layer 1 and for Gaussian
noise (`prop_gen_rank_stacksvd_singleweight_gaussian`, Track G, the same day; F18c, also
the same day, removed the hypothesis `w_i ≠ 0` that D37 had added, so the Gaussian
theorems carry no condition on the weights). The first was proved as an honest implication
from the hypothesis `hlaw` (the single-weight limit on its two-table witness at every
positive weight pair), which the discharge does not supply because the witness leaves
`EigSep` at extreme weight ratios; item F18b (2026-09-07) proves `hlaw` on the witness regime by
regime, and `prop_singleweight_suboptimality_gaussian` is unconditional. Neither is counted
among the results with no Lean statement.

The removal remark fails because removing a component with `β_ij = 0` can break
`Rank(∑_i R_i R_iᵀ) = r`, after which `W⋆` has no eigengap; in 4322 of 20000 random draws with
one `β_i = 0` plus the rank condition the gap was lost (`notes/FLAGGED.md` Q8).

**7.6 The paper's `M → ∞` remark is not formalized.** `M : ℕ` is fixed in every statement. The
remark at `main_paper.tex:407`, that `prop:stacksvd_general` survives `M → ∞` when `‖c‖₁` and
`‖θ‖₂²` converge, has no Lean form (`notes/archive/AUDIT_PACK_V4.md` S2).

**7.7 There is no bridge from `UnalignedModelR` at `rk = 1` to `UnalignedModel`.** The two are
separate structures, and no `perfRG = perfR` bridge exists, so a general-`r_i` theorem does not
visibly imply its `r_i = 1` twin. Both twins are proved separately and both are listed in the
table of section 5 (`THEOREMS.md` 1.1; `notes/FLAGGED.md`, T8 audit debt item 2).

**7.8 Two rank-`r` performance functionals are spectral surrogates, not the paper's
Frobenius norm.** `perfRG` (`RankR/General.lean:474`) and `perfRGW` (`:551`) are
`tr(Bᵀ specInvTop A r B)`; `perfStackRG` (`RankR/SubspaceG.lean:482`) is the projector sum
`∑_{k<r} ‖specProjTop(A, r) (V e_k)‖²` (corrected 2026-09-02: an earlier version of this
paragraph called it a trace through `specInvTop`). Both forms equal `‖V̂ᵀ V‖_F²` on the event
that the top-`r` eigenvalues separate from the rest and (for the trace) are positive; at a
boundary tie the spectral projector has rank above `r` and the surrogate is the larger
(`notes/archive/AUDIT_PACK_V4.md` C16, `notes/archive/audit_independent_trackB_2026-09-02.md` W1). At `r_i = 1`
both stackSVD and SVDstack have a `_frobenius_eig` twin that states the paper's own quantity;
at general `r_i` the SVDstack twin is
`prop_general_rank_unweighted_svdstack_general_frobenius_eig` and, since 2026-09-03 (F23), the
stackSVD twin is `prop_stacksvd_subspace_general_frobenius_eig_gaussian`
(`RankR/SubspaceGaussian.lean`, section 4b): the gap event `ae_topGap_stackGramG` is proved at
general `r_i` by the same route as at `r_i = 1`, and the conclusion is
`frobSq ((topEigMat ...)ᵀ * V)`, the paper's `‖V̂ᵀ V‖_F²` at the canonical top-`r` eigenframe,
under the model, the Gaussian noise and the regime only. Before that day the stackSVD twin
did not exist (`notes/archive/AUDIT_PACK_V4.md` S12, C16; both `prop_stacksvd_subspace_frobenius`
forms of `RankR/Frobenius.lean` take `UnalignedModel`, the `r_i = 1` model), and
`prop_stacksvd_subspace_general_gaussian` recovered the paper's quantity only on the gap
event. The surrogate forms stay in the tree; the Frobenius forms are the paper-literal ones.

**7.9 `w = 0` is a scalar convention, not an estimator.** `Scalars.Lw θ c 0 = 0`, while
`stackPerfW m 0 N ω` is identically `1`, because the weighted stack is the zero matrix and the
top eigenspace of the zero matrix is the whole space. Every estimator theorem therefore either
asks for a nonzero weight in its signature or asks for a `HeteroLaw` at that weight, which is
unsatisfiable at `w = 0` as soon as `2 ≤ d N`. The rank-`r` analog is `rank W ≥ r`
(`notes/archive/AUDIT_PACK_V4.md` S9, C11). The paper's own optimal stackSVD weights are all zero at
`θ ≡ 0`, where the paper's display then reads `1 → 0`; that is finding E8 and the hypothesis
`hθ : ∃ i, θ_i ≠ 0` of `thm_stacksvd_weighted_gaussian` (`notes/paper_edits.md` E8,
`notes/REVISION_LIST.md` R5).

**7.10 A Layer 1 statement alone is not evidence that anything happens.** The abstract law
structures carry no `IsProbabilityMeasure` instance, so a Layer 1 theorem also holds for the
zero measure, where `TendstoInProb` is trivially true. Every Gaussian corollary pins a
probability measure, through `GaussianNoise` or an explicit `[∀ N, IsProbabilityMeasure (μ N)]`
(`notes/archive/AUDIT_PACK_V4.md` S7, C9). Section 6 answers the matching question for the Gaussian
endpoints.

**7.11 The uniform-in-`w` bounds quantify over an event that is not shown to be measurable.**
`svdstackPerfW_uniform_bound`, `perfRW_uniform_bound` and `perfRGW_uniform_bound` have the form
`μ {∃ w, L⋆ + ε ≤ perf w} → 0`. The proofs use outer-measure monotonicity only, so the
statements are correct as written; this is recorded because a reader may expect measurability
(`notes/archive/AUDIT_PACK_V4.md` C19).

**7.12 `TopGap A hA r` is vacuous when `r ≥ p`.** No pair of indices exists there, so the
predicate holds for every symmetric `A`. That is the paper's own convention
`λ_{r̃+1} := -∞`. Do not read `hgap` as a spectral separation unless `r < p`
(`notes/archive/AUDIT_PACK_V4.md` C18, `THEOREMS.md` 6.2).

**7.13 Two hypothesis fields are stronger than anything the paper assumes, and one field is
the weighted conclusion itself.** The first two are strengthenings, not gaps, and the Gaussian
discharge proves all three (`THEOREMS.md` 7.1, `notes/archive/AUDIT_PACK_V4.md` S6):

- `SingleTableLaw.delocUniform` and `TableLawR.delocUniform` are uniform over the test
  direction; the paper's part 2 is sequential (the sequential form is the theorem
  `SingleTableLaw.deloc_seq`).
- `SingleTableLaw.topSimple` and `TableLawR.simple` assert finite-`N` almost-sure simplicity,
  which the paper does not state; they replace the paper's implicit "the top singular vector
  is well defined".
- `HeteroLaw.align` is not a strengthened assumption of the paper: it is the weighted random
  matrix conclusion (the closed form `L(w)`, one deterministic step past the BEJ1129 display
  the paper cites at line 1414, decision D21) placed behind a Layer 1 interface. So
  `thm_stacksvd_weighted_general` is a field projection, and the mathematical content of the
  weighted stackSVD chain lives in the Gaussian discharge `heteroLaw_of_gaussian`
  (`notes/archive/AUDIT_PACK_V4.md` S8, C3; external rank-one audit of 2026-09-03, section 2.7).

**7.14 Two structures deliberately omit a field.** `HeteroLaw` has no `b`-general half of
`eq:weighted_norm`, because no result consumes it; `TableLawR` has no `lamMax` field, because
no statement of Section 7 reads it (`THEOREMS.md` 2.2, 2.3).

**7.15 The MLE result: the identity, and since 2026-09-03 the derivation too.**
`MLE.mleLogLik` writes the paper's marginal log-likelihood as a formula, with `c_i` where the
paper writes `n_i/d`; `mleLogLik_eq` and `mleLogLik_max_iff_mem_topSpace` prove that the
written objective is maximized exactly on the top eigenspace of the weighted stack Gram
matrix (`THEOREMS.md` 5.5, `notes/archive/AUDIT_PACK_V4.md` C15). Item L5 of the external statement
audit asked for the Gaussian marginalization over `u_i` that produces the formula.
`thm_wstacksvd_mle_marginal` (`MLEMarginal/Main.lean`) now proves it: the random-effects
model (`u_i ~ N(0, I/n_i)`, Gaussian noise, independent tables) is the measure
`reJointLaw`; its Lebesgue density is the explicit `reDensity`
(`∏_i ∏_k gaussDensity d (Σ_i(v)) (X_i)_k`); `log reDensity` at the observed tables equals
`mleLogLik` at `c_i = n_i/d` up to `(∑_i n_i) d/2 · log(2π)`; and its unit maximizers are
the same top eigenspace. What remains a modeling choice: `reJointLaw` is a separate object,
and nothing is claimed about the law of the fixed-effects model `m`
(`notes/archive/L5_mle_marginal.md`). Finding E5.3, that the paper should name the
marginalization as a standard step, stands as a presentation remark.

**7.16 The degenerate case `A_β = I` is excluded.** When exactly one `β_i > 0` the Python and
the R reference code disagree (`β₁²` against `0`), Monte Carlo gives about `0.33`, and the
paper leaves the case open. Lean excludes it through `hthr` (`notes/FLAGGED.md` Q2,
`THEOREMS.md` 7.3).

**7.17 Everything is convergence in probability, and this is not a gap.** The paper claims the
same mode (`main_paper.tex:234`, `:159`, `:175`); section 3.8 gives the evidence. Nothing is
stated almost surely in `N`; the two almost-sure facts (`SingleTableLaw.topSimple`,
`TableLawR.simple`) are per-`N` statements, not limits (`THEOREMS.md` 7.1,
`notes/archive/AUDIT_PACK_V4.md` S3).

**7.18 Two paper labels an auditor will grep and not find.** `prop:maxrow`
(`main_paper.tex:2628`) and `cor:domination` (`:2679`) are commented out in the source of the
AoS resubmission, together with the `\maxrow` estimator. They are not claims of the paper and
nothing formalizes them (`THEOREMS.md` 7.1, `notes/README.md`).

**7.19 One paper lemma has no section of its own.** `lem:noise_projection_concentration`
(`main_paper.tex:2496`) is proved, not assumed, as `SpikedModel.noise_projection_tendsto`
(`ThetaEst.lean`), and only `thm:theta_est` reads it (`THEOREMS.md` 7.1).

**7.20 Well-formedness hypotheses that remain.** These make an object exist or an index be in
range. An auditor should read each one and check that it does no hidden work
(`THEOREMS.md` 7.1):

| Hypothesis | Where | What it does |
|---|---|---|
| `[∀ N, IsProbabilityMeasure (μ N)]` | every Gaussian theorem | each `μ N` is a probability measure; derivable where `JointGaussianNoise` is present, and several theorems do without it |
| `[NeZero M]`, `hM : 0 < M` | the rank-`r` results, and the 26 rank-one declarations quoted in `THEOREMS.md` that sit inside a section or namespace with `variable [NeZero M]` in force | at least one table, so the stack has rows (`MultiTableModel.stack` reads `(m.tbl 0).v`, `StackSVD.lean:157`). In the rank-one files the instance is a section-level `variable [NeZero M]` (13 files, listed in `THEOREMS.md` section 0 item 3), which Lean adds to every declaration that mentions `M`; until 2026-09-02 the quoted source text of 26 declarations did not show it while their elaborated types carried `[inst_1 : NeZero M]` (external statement audit finding 11, which counted 21; `#check` output in `notes/archive/audit_external_statement_nezero_check_2026-09-02.txt`). Since item F22 (2026-09-02) each of the 26 writes `[NeZero M]` as its first binder, with `omit [NeZero M] in` in front so the section instance is not duplicated; `notes/archive/F22_nezero.md` lists them |
| `[NeZero S.card]`, `S.Nonempty` | `cor.2`, `prop:binarystacksvd_inadmissable` | the kept subset is nonempty |
| `hr : 0 < r` | the rank-`r` theorems that read a top-`r` eigenframe | the shared subspace is nontrivial; at `r = 0` both sides are `0`. The 2026-09-02 cleanup removed it from 8 declarations that never read it and kept it at the other 25 |
| `hrr : r ≤ rtot rk`, `hrM : r ≤ M` | sections 6.2, 6.4 of `THEOREMS.md` | the paper's `r ≤ r̃` |
| `hc : 0 < c i` | almost everywhere | `eq:RMT_limit`'s `c_i ∈ (0, ∞)`; `prop:stacksvd_subspace` needs only `0 < ∑ i, c i` |
| `hβdef`, `hθ`, `hR` | many | naming hypotheses; they tie a free variable to the model and prove nothing |

**7.21 One document in the packet's evidence trail is a snapshot and is stale in two places.**
`notes/archive/AUDIT_PACK_V4.md` is dated 2026-09-02 03:10. Its scope item S4 (rank-`r` `TableLawR`
assumed at `r_i ≥ 2`) and S10 (`thm:rank_r_stacksvd` has no Lean statement) were closed later
the same day by Track C and Track E (`notes/PROGRESS.md`, "Where we are" 11:00 and 12:30;
`TECHNICAL.md`). Items S1, S2, S6 to S9, S11, S12 and the caveats C1 to C21 are still current and
are the source of several rows above.

**7.22 One packaging defect.** Every Lean file header says "Released under Apache 2.0 license
as described in the file LICENSE", and the repository tracks a `LICENSE` file only under
`lean/StackedSVD/Vendor/COLT83/` (`git ls-files`, run this session). The project's own
LICENSE file is absent.

---

## 8. Findings about the paper, E1 to E12

The paper is **unchanged**. Claude never edits it (`CLAUDE.md` hard rule 1). Every item below
is a proposal for the authors, and the status of all twelve (E12 added 2026-09-07, with F18b) is `proposed`; the
user has applied none of them (`notes/paper_edits.md`, opening table; `THEOREMS.md` 7.2). The user's
decision of 2026-09-05: the paper is out for review, so the list stands as the clarifications
the authors will make in a later revision, and the external audit reads the paper as it is,
with this list next to it. Each names
a Lean anchor, so the claim is checkable. `notes/paper_edits.md` opens with the same twelve items in
one line each and then gives, for each finding, the exact old and new LaTeX of every proposed
edit (29 replacements, 24 on 2026-09-05, 4 on 2026-09-06 and 1 on 2026-09-07, subsections `Exact replacement Ek.j`; since the
snapshot of 2026-09-09 each carries a `Status:` line, `open`, `applied ..., verbatim.` or `applied ..., with changes.`, the
last with an `As applied` block that quotes the text the user wrote instead); `THEOREMS.md`
quotes each amended statement environment under `Paper, as amended`, and the sixth gate
`scripts/check_paper_edits.py` checks both against the snapshot (section 9). The kind of evidence differs by item: E1, E2, E3,
E4, E5, E8 and E9 are backed by Lean theorems or by direct reading of the definitions; E6 is
a closed-form check; E7 is a closed-form outlier computation plus Monte Carlo (7.4); E10 is
a bibliographic check; E11 is a closed-form computation plus Monte Carlo, now also backed by
the Lean anchor `EigSep.sorted` (`THEOREMS.md` 6.8, added 2026-09-05).

**E1. `thm:svdstack_weighted` (line 502): optimality holds over all weights, with no eigengap
condition.** The paper says "an optimal weighting of SVDstack is `w⋆`" and its proof starts
with "Provided that `W A_β W` has a unique largest eigenvalue", so it compares `w⋆` only with
weights whose limit matrix has an eigengap; `w = (1, 1, 0, …, 0)` with `β_1 = β_2 = 0` (top
eigenvalue `1` of `diag(1, 1, 0, …)` with multiplicity `2`) and any data-dependent weight are outside
that class. The evidence is a Rayleigh-quotient bound: every unit `x` in the top right singular
subspace of `W Ṽ` is `Ṽᵀ b`, so `⟨x, v⟩² = (bᵀg)²/(bᵀGb) ≤ gᵀG⁻¹g`, and the right side does not
depend on `w`. The proposal is to replace the optimality sentence by a `limsup` bound valid for
every weighting, deterministic or data-dependent, and to delete or reduce the footnote at line
888. Anchor: `svdstackPerfW_le_rowBound` (finite `N`, almost sure), `rowBound_tendsto`,
`svdstackPerfW_uniform_bound`, and `thm_svdstack_weighted_gaussian_opt_full`
(`SVDStack/Rayleigh.lean`).

**E2. `thm:gen_rank_weight_svdstak` (line 893): drop `β_ij > 0` and the footnote.** The
hypothesis "`β_ij > 0` for all `i, j`" excludes every table that does not see every component,
which is the setting the unaligned model is built for, and the remark after the theorem (a
component with `β_ij = 0` "can be removed") changes `r̃` and can make `rank B_R < r`. The
evidence is the same row-space bound at rank `r`, `‖V̂(W)ᵀV‖_F² ≤ tr(gᵀG⁻¹g)`, plus the identity
`tr(B_Rᵀ A⁻¹ B_R) = L⋆` at every rank of `B_R` (numeric agreement `1e-9`), plus attainment at
`W⋆ = D^{-1/2}`. A prior instance with `rank B_R = 1` at `d = 400, 800, 1600` gave means
`0.420, 0.402, 0.395` against `L⋆ = 0.400`. Step 2 of the printed proof (`:2062`) reads
"`diag(β)` is invertible by assumption", the dropped hypothesis, so E2.4 (added 2026-09-06)
rewrites it: the bound at every `N`, the trace identity at every rank of `B_R`, and the
attainment at `W⋆` with or without the eigengap; the `rank B_R < r` half of that argument has
no Lean twin (the Lean attainment without a rank hypothesis is the trace form,
`thm_gen_rank_weight_svdstak_general_r_norank`). The proposal is to remove the positivity
hypothesis and to replace the footnote by the uniform bound with the single requirement
`rank W ≥ r`. Anchors: `RankR/WeightedUpper.lean`, `RankR/WeightedUpperG.lean`
(`perfRGW_le_rowBoundRG`, `thm_gen_rank_weight_svdstak_general_r_full`).

**E3. `eq:stacksvd_gammak` (line ~2331): define `γ_j` below the threshold.** The paper defines
`γ_j` as the unique root in `(0,1)` of `Σ_i θ_ij⁴(1−x)/(c_i + xθ_ij²) = 1`. A root exists if
and only if `Σ_i θ_ij⁴/c_i > 1`, so below that threshold `γ_j` is undefined and
`thm:rank_r_stacksvd` has no meaning for component `j`. The proposal is to add "and `γ_j = 0`
otherwise", which matches the rank-one weighted result. Anchor: `Scalars.gammaR`, total by
construction (`Scalars.stackSVDLimitW` is `0` below the threshold).

**E4. `thm:rank_r_svdstack` (line ~2400): it is not a specialization.** The corollary is
labeled a specialization of `thm:gen_rank_weight_svdstak`, but that theorem as stated requires
`β_ij > 0` and the exactly aligned model has `β_ij = 0` at every subcritical pair; the
component clause also needs the `j`-th eigenvalue of the block-diagonal limit to be simple. The
evidence is that the component clause needed its own proof in Lean. The proposal is to derive
the aggregate clause from E2 and to state the component clause under its own hypotheses. One
sub-item of E4 should itself be revised: item 2 asks for `S_j > 0`, and the clause holds at
`S_j = 0` with limit `0` (necessity scan, seed 20260904: `(v_1ᵀ v̂_1)²` falls from `0.0277` at
`d = 200` to `0.0011` at `d = 1600`). Anchors: `thm_rank_r_svdstack_component_gaussian`,
`thm_rank_r_svdstack_aggregate_gaussian` (`notes/paper_edits.md` E4, `notes/FLAGGED.md` item 14).

**E5. Minor items, four of them.** `thm:simple_thm1` is valid for every `M ≥ 1` and the two
branches agree at the threshold (anchor `thm_simple_thm1_svdstack_gaussian_full`); the strict
clauses of `thm:stacksvd_binary_optimal_svd_stack` and `prop:dominance` need no change, and are
listed only so the Lean strict wrappers have a paper location; the MLE appendix should name the
Gaussian marginalization as a standard cited fact, so that the formal and the informal halves
are separate (the step itself is proved in Lean since 2026-09-03, 7.15); the two remarks need
no change.

**E6. `remark:svd_outperform_stack` (line 603): the binary example sits on the threshold.** The
paragraph takes `θ₃ = c₃^{1/4}`, so `θ₃⁴ = c₃`, and `cor.2` defines the binary set with a
**strict** inequality, so table 3 is discarded. Binary-weighted stackSVD is then the stack of
tables 1 and 2 with limit `(8² − 2)/(8² + 8) = 62/72 = 0.8611`, above SVDstack's `6/7 = 0.8571`,
the reverse of the claim. The evidence is a deterministic closed-form check
(`scripts/numeric/check_remark_examples.py`, PASS). The proposal is `θ₃ = (c₃ + 1)^{1/4}`, which gives
binary stackSVD `0.1413` against SVDstack `0.8570` at `c₃ = 10⁴`. Anchor:
`remark_svd_outperform_stack_three_binary`, which states the true comparison.

**E7. `thm:rank_r_stacksvd` (line ~2306) and `alg:rank_r_stacksvd` (line 2382): the index `ℓ_j`
is wrong at unordered `θ`.** The algorithm sets `ℓ_j` to the rank of `θ̃_jj` among `{θ̃_jk}_k`,
but under block heteroscedastic noise the outlier of spike `k` depends on the whole profile
`(w_ij² θ_ik²)_i`, not on its sum. On the instance `θ = [[3.61993204, 0.40121508], [0.09540074,
1.52018581]]` with `c = (1, 1)`, the paper's index is 1 while the top singular vector carries
spike 0, and Monte Carlo from `d = 200` to `1600` sends the overlap at index 0 to `0.57` and at
the paper's index 1 to `0`. The proposal is to define `ℓ_j` by the rank of the outlier `ρ_jk`,
or to state the same-ordering hypothesis, under which the two definitions agree. Anchors:
`Scalars.ellR`, `Scalars.ellSup`, `Scalars.ellSup_eq_ellR` (`RankR/Het/Scalars.lean`).

**E8. `thm:stacksvd_weighted` (line 463) and `prop:dominance` (line 634) are false at
`θ ≡ 0`.** The stated optimal weights `w_i⋆ ∝ θ_i/√(θ_i² + c_i)` are all zero there, the
weighted stack is the zero matrix, its top eigenspace is the whole space, and `(vᵀ v̂)²` is
`1` at every `N`, while the theorem's own "otherwise `γ⋆ = 0`" branch claims `0`. The proposal
is the hypothesis `∑_i θ_i² > 0`. Anchor: `hθ : ∃ i, θ_i ≠ 0` in `thm_stacksvd_weighted_gaussian`
and `prop_dominance_gaussian` (`RMT/Het/Sup.lean`), which is why those signatures carry it
(`notes/REVISION_LIST.md` R5).

**E9. Eight errata and degenerate cases** (`main_paper.tex:376, 442, 1331, 1612, 842, 909, 1345, 212`): `M ≥ 2`
is missing where `M = 1` makes a statement empty; a maximum runs over a set that contains the
empty subset; a weight condition describes an empty set; an inequality has its sign reversed;
`V̂_stacksvd` is not unique at a tie; two counts `2^M` include the empty weighting; a caption
cites `thm:stacksvd_weighted` for `S`, which `thm:svdstack_weighted` defines. Detail and anchors: `notes/REVISION_LIST.md` R1 to R4 and
R8.

**E10. The heteroscedastic citation points to a white-noise paper** (`main_paper.tex:403,
1372`): the paper describes Theorem 2.3 of `10.3150/19-BEJ1129` as a result on "signal-plus-noise
matrices with heteroscedastic noise". That reference (Ding, Bernoulli 26(1), 2020,
arXiv:1702.06975) assumes white noise (Assumption 1.1, i.i.d. entries of variance `1/N`) and
does not contain the words "heteroscedastic" or "variance profile". Kind of evidence: a
bibliographic check of the ar5iv rendering on 2026-09-02 (identity, Assumption 1.1, Theorem
2.3, full-text word search); the published PDF was not read. The Lean tree does not import the
heteroscedastic input: `heteroLaw_of_gaussian` (`RMT/Het/Sup.lean`) proves it for Gaussian
noise (section 5 of this document, `THEOREMS.md` section 4).

**E11. The distinctness clause of `assum:gen_rank_stacksvd_eig_sep` is needed for the
value** (`main_paper.tex:2104`, `:2112`; added 2026-09-05): at a tied root `γ_ℓ = γ_m` the
eigenvectors `z_ℓ, z_m` of `prop:gen_rank_stacksvd_singleweight` are defined only up to a
rotation inside a plane, and the sum of reciprocals in the performance formula is not
invariant under it. The proof's own construction fixes the choice: the `z_ℓ` of a tied plane
must be orthogonal for the form `K(γ_ℓ) = ∑_i (w_i² / (w_i² - γ_ℓ)²) R_i Θ_i² R_iᵀ`. A remark,
not an error, since the assumption excludes ties. Kind of evidence: closed form plus Monte
Carlo (`scripts/numeric/check_singleweight.py`, seed 20260905): on `M = 2`, `r = 2`, `θ² = (3, 15)`,
`w = (1, 0.5)`, `c = (0.05, 0.05)`, the `K`-eigenbasis gives 1.677750 against Monte Carlo
1.677916 ± 0.004956 at `d = 1600`, and a 45 degree rotation gives 1.657037 (off by 4.2 to 6.2
standard errors). Lean anchor: `EigSep.sorted : StrictAnti γ`
(`RankR/SingleWeight/Scalars.lean`, `THEOREMS.md` 6.8), the field that turns this distinctness
clause into a strict order and rules the tie out by hypothesis, added 2026-09-05.

---

## 9. How to check it yourself

Order matters: the axiom gate reads the compiled oleans, so build first
(`scripts/README.md`, the stale-oleans caveat).

```sh
cd lean && lake exe cache get && lake build   # rechecks every proof
cd .. && scripts/check_sorries.sh
CHECK_AXIOMS_BUILD=0 scripts/check_axioms.sh
python3 scripts/check_root_imports.py
python3 scripts/check_theorems_sigs.py docs/THEOREMS.md
scripts/check_layering.sh
python3 scripts/check_paper_edits.py
```

(On the development server a capped wrapper, `lake-build-capped.sh 6`, stands in for bare `lake build`.)

| Step | Runtime | Expected result |
|---|---|---|
| root `lake build` | 27 minutes on 6 cores for the full compile of 2026-09-07 (161 files, the slowest 83 s); 40 minutes on 6 cores for a full compile (coordinator, 2026-09-02); `notes/PROGRESS.md` records 35 minutes on 6 cores for the 12:21 build of 2026-09-02, 24 minutes for the 16:06 build with 107 jobs rebuilt, 25 minutes for the 01:01 build of 2026-09-03 with 7 jobs rebuilt (`Prob/GaussianDensity` alone took 1252 s), and 1 to 5 minutes when every job replays from cache | 2026-09-10 15:58:20 to 16:01:47 EDT, `Build completed successfully (8978 jobs)`, exit 0, 0 errors, 6 built (`RMT.General.Edge.CountBound`, `Edge.Excess`, `Edge.Markov`, `Edge.Sup`, `General.Layer1`, the root `StackedSVD`), the rest replayed (built by the root build of 15:53:55 to 15:57:59 EDT on the same sources up to one docstring line, exit 0, 8978 jobs, which compiled `Edge.Dyck` in 42 s, `Edge.Code` in 58 s and `Edge.CountBound` in 32 s, and by the integration unit's builds of 15:40 to 15:46 EDT), on the committed tree `416ab9d` (worktree clean); the transcript of Appendix A is the audit pack header run of 16:18:21 to 16:19:22 EDT on `416ab9d`, 0 rebuilt, 5 replayed, 61 s, worktree clean. Before that: 2026-09-10 12:37:28 to 12:39:20 EDT, `Build completed successfully (8975 jobs)`, exit 0, 0 errors, 3 built (`RMT.General.Edge.Sup`, `General.Layer1`, the root `StackedSVD`), the rest replayed (built already by the targeted builds of `Edge/Count.lean`, `Edge/Excess.lean` and `Edge/Markov.lean` at 12:34 to 12:36 EDT, `Build completed successfully (8782 jobs)`, exit 0, and by the root build of 11:04 EDT on the earlier Stage 3 files), on the working tree of the Stage 3 commit that follows the stub commit `aa27519` on `stage3` (not `aa27519`'s own committed tree: the eleven `RMT/General/Edge/` files, `General/Layer1.lean` and the root import list were uncommitted at build time); the transcript of Appendix A is the audit pack header run of 13:04:53 to 13:05:47 EDT (UTC) on `aa27519`, in that same uncommitted state, 0 rebuilt, 6 replayed, 54 s; its worktree line lists 22 files modified or untracked at header time, 12 of them Lean, so the Lean tree of the run is the Stage 3 working tree, not `aa27519` exactly. Before that: 2026-09-10 09:02:26 to 09:17:10 EDT, `Build completed successfully (8964 jobs)`, exit 0, 0 errors, 112 built (the closures of the 10 files whose docstrings cite moved documents, `StackSVDWeighted.lean` among them; 21 s to 90 s each), the rest replayed, on the Lean tree of the rename commit `baf31a4` (the lake folder `StackedSVD/` renamed to `lean/`; comment text only, no module name and no statement changed; the type dump of every constant, 7987 rows, byte-identical to that of `6633055`); the transcript of Appendix A is the audit pack header run of 09:33:16 to 09:34:15 EDT on `66be99d` (the documents commit on top of `baf31a4`, no Lean file), 0 rebuilt, 5 replayed, 59 s; its worktree line lists the 2 Markdown files of the docs commit that follows. Before that: 2026-09-10 07:32:33 to 07:56:37 EDT, `Build completed successfully (8964 jobs)`, exit 0, 0 errors, 197 built (all 196 files of ours and the root, since the edited files reach nearly every module; 21 s to 101 s each, `ResolvDeriv` 101 s), the 20 vendored files replayed, on the Lean tree of the docstring commits `848f7d0` and `99afeab` (89 stale `notes/<name>.md` paths in 126 files now `notes/archive/<name>.md`, the 15 citations of session scratch scripts rewritten; comment text only, the type dump of the 7987 constants byte-identical to that of `da974e9`); the root-only rebuild of 08:19:33 to 08:20:12 EDT on `6633055` (two comment lines of the root file, the last stale note paths) built 1, the root, 22 s, its olean byte-identical; the transcript of Appendix A is the audit pack header run of 08:23:32 to 08:24:23 EDT on `6633055`, 0 rebuilt, 5 replayed, 51 s; its worktree line lists the 5 Markdown files of the docs commit that follows, no Lean file. Before that: 2026-09-09 23:02:21 to 23:11:20 EDT, `Build completed successfully (8964 jobs)`, exit 0, 0 errors, 19 built (the 14 edited files of `RMT/General/`, their 3 readers there, `General/Layer1`, the root; 21 s to 55 s each, `Companion` 201 s under NFS load), the rest replayed, on the Lean tree of the cleanup commit `da974e9` (F36 to F38: 50 duplicate private helpers deleted, 38 made public, no theorem statement changed; the coordinator's replays of 23:26:37 to 23:27:38 EDT, 1 built, `RMT/General/Iso`, olean unchanged, and 23:28:22 to 23:28:32 EDT, 0 built; the transcript of Appendix A was the audit pack header run of 23:41:32 to 23:42:37 EDT on `b8771cb`, 0 rebuilt, 5 replayed, 65 s, its worktree line the 2 Markdown files of the docs commit that follows). Before that: 2026-09-09 21:23:10 to 21:29:03 EDT, `Build completed successfully (8964 jobs)`, exit 0, 0 errors, 24 built (the closure of `RMT/General/Defs.lean`: the 22 modules of `RMT/General/` and `Prob/Chebyshev` that import it, `General/Layer1`, the root; 19 s to 50 s each), the rest replayed, on the Lean tree of the F35 commit `d584cc6` (`ResolventFormsC` moved from the retired `RMT/General/Statements.lean` to `RMT/General/Defs.lean`, no statement changed; the coordinator's replay of 21:33:55 to 21:34:09 EDT built 0 jobs; the transcript of Appendix A is the audit pack header run of 21:51:14 to 21:52:10 EDT on that tree, 0 rebuilt, 5 replayed, 56 s; its worktree line lists the 13 Markdown files of the docs commit that follows, no Lean file). Before that: 2026-09-09 20:26:16 to 20:30:00 EDT, `Build completed successfully (8965 jobs)`, exit 0, 0 errors, 10 built (`RMT/General/Statements.lean` 35 s, `FormsSup` 21 s, `DelocLimits` 32 s, `FormsGeneral` 27 s, `DelocAlign` 35 s, `DelocUniform` 45 s, `DelocUniformSup` 23 s, `Sup` 24 s, `General/Layer1` 24 s, the root 23 s), the rest replayed, on the Lean tree of the Stage 1 close commit `304bffc` (units G4, F1a, F1b, F2, F3, D1a, D1b, D1c, G7: 13 new files under `RMT/General/`, the 12 target statements all discharged, no `sorry` left; the 11 files of wave 5 were compiled by the `lake-build-capped.sh 12` run of 19:51:26 to 19:54:20 EDT, 8963 jobs, exit 0, 8 built, on the wave 5 tree with the one stub, not committed in that state; the transcript of Appendix A is the audit pack header run of 20:43:58 to 20:45:07 EDT, 0 rebuilt, 5 replayed, 69 s on the committed tree). Before that: 2026-09-08 22:53:23 to 23:07:37 EDT, `Build completed successfully (8939 jobs)`, exit 0, 0 errors, 109 built (the closures of `StackSVD.lean`, whose header comment changed, `StackSVD/Main.lean` and `RMT/Het/Simplicity.lean`; 18 s to 70 s each), 5 replayed (the vendored COLT83 files, the only jobs with a stored log), on the Lean tree of the F32 commit `9c33ae6` (the two duplicate private helpers dropped; the transcript of Appendix A is the audit pack header run of 23:35:44 to 23:36:36 EDT on that committed tree, 0 rebuilt, 5 replayed, 52 s). Before that: 2026-09-08 21:37:13 to 21:39:23 EDT, `Build completed successfully (8939 jobs)`, exit 0, 0 errors, 5 built (`RankR/SubspaceGaussian.lean` 23 s, `RankR/SingleWeight/Tie.lean` 22 s, `RankR/SingleWeight/Suboptimality.lean` 24 s, `Main.lean` 24 s, the root 20 s), the rest replayed, on the Lean tree of the F28 commit `16c100c` (the import boundary of the Gaussian RMT core; header run of 21:57:30 to 21:58:24 EDT on that committed tree, 0 rebuilt, 5 replayed, 54 s). Before that: 2026-09-08 15:54:40 to 16:32:30 EDT, `Build completed successfully (8933 jobs)`, exit 0, 0 errors, 6 built (`Prob/NoiseLaw.lean` 1871 s and `Prob/LinFormMoments.lean` 1879 s, both in the NFS import phase, `Prob/NoiseMoments.lean` 28 s, `ThetaEst.lean` 38 s, `Main.lean` 156 s, the root 23 s), 5 replayed (the vendored COLT83 files, the only jobs with a stored log), on the Lean tree of the Stage 0 commit `9fdcbfa` (the transcript of Appendix A is the audit pack header run of 16:50:06 to 16:55:57 EDT on that committed tree, 0 rebuilt, 5 replayed, 351 s under NFS load). Before that: 2026-09-08 08:14:24 to 08:16:49 EDT, `Build completed successfully (8930 jobs)`, exit 0, 0 errors, 3 rebuilt (`ThetaEst.lean` 39 s, `Main.lean` 34 s, the root 32 s), 5 replayed, on the Lean tree of commit `8736b2a` (F31; header run of 08:32:05 to 08:32:56 EDT on that committed tree, 0 rebuilt, 5 replayed). Before that: 2026-09-07 21:55:38 to 22:16:22 EDT, `Build completed successfully (8930 jobs)`, exit 0, 3 rebuilt (`Secular.lean` 1125 s under NFS load), on the tree of `18452a2`; 2026-09-07 18:14:06 to 18:34:34 EDT, `Build completed successfully (8929 jobs)`, exit 0, 0 errors, 114 rebuilt (the closure of the last two `omit` edits of the linter campaign, 18 s to 68 s each), 5 replayed (the vendored COLT83 files, the only jobs with a stored log), on the Lean tree of commit `926701c`; the other edited modules were compiled by the `lake-build-capped.sh 6` runs of 17:24:18 to 17:51:15 EDT (8929 jobs, 161 rebuilt, 27 min, the slowest file `RMT/ResolvDeriv` 83 s) and 17:52:22 to 18:10:11 EDT (8929 jobs, 117 rebuilt), both exit 0 (the transcript of Appendix A is the audit pack header run of 21:06:38 to 21:08:08 EDT on the committed tree `a32907a`, 0 rebuilt, 5 replayed, 90 s). Before that: 2026-09-05 19:23:32 to 19:24:29 EDT, `Build completed successfully (8922 jobs)`, exit 0, 0 errors, 1 rebuilt (the root, 40 s), 73 replayed, on the Lean tree of commit `85e6a0e` (F18c); the two changed modules `RankR/SingleWeight/Het/Sub` (41 s) and `Het/Sup` (34 s) were compiled by the `lake-build-capped.sh 4` run of 19:21:28 to 19:23:07 EDT (8875 jobs, exit 0) and replay from the cache since; the six other `RankR/SingleWeight/Het/` modules were compiled by their own `lake-build-capped.sh 4` runs of 15:13 to 16:29 EDT (`Forms` 8844 jobs, `Outliers` 8862, `Align` 8867, all exit 0) and the root build of `a2cbbb9` (16:43:55 to 16:45:50 EDT, 8921 jobs, 9 rebuilt, exit 0); the 146 modules before them were compiled by the builds of `947b453` (13:26 to 13:41 EDT), `482309c` (12:20 EDT), `697da9c` (09:57 to 10:37 EDT), `a5decd4` (2026-09-03 15:33 EDT), `04fd96b` (09:02 EDT) and `e2c76c5` (01:01 to 01:26 EDT), all exit 0 (the transcript of Appendix A is the audit pack header run of 19:36:07 to 19:36:57 EDT on the committed tree, 0 rebuilt, 73 replayed, 50 s) |
| `scripts/check_sorries.sh` | 1.5 s (measured this session) | `check_sorries: OK - 0 sorry sites in 0 declarations, 0 rows, all matched.`, exit 0 (1 sorry site in 1 declaration, 1 row at `0c52d54`, before the X8 count discharged it: `cellCard_mul_le_pos`, `RMT/General/Edge/Count.lean`, the Furedi-Komlos shape count, tracked in `docs/SORRIES.md`, caveat 8 of `README.md`; 0 sorry sites in 0 declarations at `aa27519`/`baf31a4`, before Stage 3; the 12 target statements of the non-Gaussian Stage 1 held a tracked `sorry` each from the morning of 2026-09-09 to the evening of the same day, the last of them `delocUniform_of_general` until 20:25 EDT; no build of record was taken in that state) |
| `CHECK_AXIOMS_BUILD=0 scripts/check_axioms.sh` | a few minutes over NFS (coordinator, 2026-09-02). measured on the run of record: under 1 minute, from the end of the build to the end of the gate. `notes/archive/audit_docs_2026-09-02.md` check 8 reports 15 to 25 minutes for a cold NFS run under load | `check_axioms: 5631 declarations audited (5631 3297).` (5327, 2789 at `0c52d54`, before the X8 count; 5121, 2600 at `aa27519`/`baf31a4`, before Stage 3; 5171, 2603 at `d584cc6`, before the F36 to F38 dedup deleted 50 private helpers; 5171, 2601 at `304bffc`; 5160, 2598 on the wave 5 tree of 2026-09-09 with the one stub, not committed; 4687, 2346 at `9c33ae6`; 4689, 2348 at `16c100c`; 4688, 2341 at `9fdcbfa`; 4622, 2329 at `8736b2a`; 4618, 2324 at `18452a2`; 4576, 2323 at `926701c`; 4519, 2249 at `9c82538`; 4500, 2243 at `a2cbbb9`; 4361, 2139 at `947b453`; 4261, 2084 at `697da9c`) then `check_axioms: OK - no axiom outside {propext, Classical.choice, Quot.sound} and no untracked sorry.`, exit 0 |
| `python3 scripts/check_root_imports.py` | 0.2 s (measured this session) | `check_root_imports: 230 modules (210 ours, 20 vendored); 229 reachable from lean/StackedSVD.lean.` (227, 207, 226 at `0c52d54`, before the X8 count; 216, 196, 215 at `aa27519`/`baf31a4`, before Stage 3; 217, 197, 216 at `304bffc`; 191, 171, 190 at `9c33ae6`; 185, 165, 184 at `9fdcbfa`; 182, 162, 181 at `8736b2a`; 181, 161, 180 at `926701c`; 174, 154, 173 at `9c82538`), one vendored unreachable module named, `OK`, exit 0 |
| `python3 scripts/check_theorems_sigs.py docs/THEOREMS.md` | 1.1 s (measured this session) | `signatures 148, verbatim 136, no verbatim match 12`, then the four-line legend, exit 0 (140, 128, 12 at `aa27519`/`baf31a4`, before Stage 3; 133, 121, 12 at `304bffc`, before the seven blocks of the non-Gaussian theorems, F40; 128, 116, 12 at `8736b2a`; 126, 114, 12 at `18452a2`; 125, 113, 12 at `9c82538`; 121, 109, 12 at `947b453`; 97, 85, 12 at `697da9c`; 96, 84, 12 before the F23 block of 2026-09-03; 92, 80, 12 before `Existence.lean`; 71, 59, 12 before the 16 blocks added to `THEOREMS.md` on 2026-09-02; 87, 75, 12 before the 4 blocks of the L2 and L3 additions of the same day; 91, 79, 12 before the L5 block of 2026-09-03) |
| `python3 scripts/check_paper_edits.py` | 0.1 s (measured 2026-09-05); reads no Lean | one line per hunk (`E5.1   main_paper.tex:338  old 1 lines -> new 2 lines` first, 29 lines in paper order), then `check_paper_edits: PASS (29 replacements: 2 applied verbatim, 11 applied with changes, 16 open; 6 amended blocks)`, exit 0. An open replacement must occur once in the snapshot; an applied one must be gone, and its new text (or its `As applied` block) must occur once; the amended text is the snapshot with the open replacements applied (before 2026-09-09: `OK - 29 replacements, each old text found once, no overlaps, new texts balanced; 6 amended blocks in docs/THEOREMS.md match.`; 28 replacements at `9c82538`) |
| `scripts/check_layering.sh` | 30 s (measured 2026-09-03; the driver imports six modules and walks three closures) | `SELFTEST 5 consumers recognized, 3 endpoints clear`, 3 `ENDPOINT` lines, 98 `USES` lines (97 before the F32 dedup of 2026-09-08, which routed `heteroLaw_of_gaussian` through the public `MultiTableModel.sum_orderEmb`; 98 on the run of 2026-09-09 21:49 EDT), no `FORBIDDEN` or `CONSUMER` line, `LAYERAUDIT_END 3 0`, then `check_layering: OK - 3 endpoints audited, 0 forbidden hits.`, exit 0 |
| `scripts/check_kernel.sh` (gate 7, added 2026-09-07) | 686 s with 3 workers on the root-disk olean copy (measured 2026-09-10 16:04:23 to 16:15:50 EDT on the committed tree `416ab9d`; the 10-minute memory poll read 23.5 GB total RSS for the user at 16:09 EDT, during the run, and 6.1 GB at 16:19 EDT, after it; rerun on the same tree at 18:14:24 to 18:25:55 EDT after the script gained its portable mode for machines without the nprocs shim, `scripts/KernelReplay.lean`: 690 s, the same verdict, 20.5 GB total RSS at the one poll; the portable mode itself, measured the same evening in a fresh clone of the public copy, replayed its 230 oleans in 737 s with 3 workers and a peak of 20 GB). Before that: 669 s with 3 workers (measured 2026-09-10 12:41:28 to 12:52:38 EDT on the Stage 3 working tree `0c52d54`, RSS not polled this run; 680 s and RSS 13.3 GB for the largest worker at the one poll on 2026-09-10 09:19:10 to 09:30:31 EDT on the tree of `baf31a4`; 692 s and RSS 13.9 to 16.5 GB on 2026-09-10 08:06:36 to 08:18:09 EDT; 612 s and RSS 14.9 to 19.3 GB on 2026-09-09 23:30:40 to 23:40:53 EDT; 644 s and RSS 12.7 to 16.7 GB on 2026-09-09 21:36:37 to 21:47:22 EDT; 627 s and RSS 15.2 GB on 2026-09-09 20:32:27 to 20:42:55 EDT; 623 s and RSS 15.6 to 18.3 GB on 2026-09-09 20:04:41 to 20:15:07 EDT; 564 s on 2026-09-08 23:14 to 23:23:33 EDT with RSS 12.0 to 24.8 GB at the polls; 551 s and RSS 16 to 24 GB on 2026-09-08 21:42:59 to 21:52:11 EDT; 628 s on 2026-09-08 16:35:02 to 16:45:30 EDT with RSS 14 GB at the one poll, 20 GB in total; 492 s and RSS 14 to 17 GB at three polls on 2026-09-08 08:20:53 to 08:29:05 EDT; 507 s and 19.7 GB on 2026-09-07); over NFS the import of the Mathlib closure alone takes 5 to 10 minutes per module | one `replaying <module>` line per module, then `check_kernel: OK - 231 modules replayed through the kernel, 0 problems; every one of the 231 oleans of the package (686 s, 3 workers).`, exit 0 (228 modules, 669 s at `0c52d54`, before the X8 count; 217 modules, 644 s at `aa27519`/`baf31a4`, before Stage 3; 218 modules, 627 s at `304bffc`; 216 modules, 623 s on the wave 5 tree of 2026-09-09 with the one stub, not committed; 192 modules, 564 s at `9c33ae6`; 192 modules, 551 s at `16c100c`; 186 modules, 628 s at `9fdcbfa`; 183 modules, 492 s at `8736b2a`). The gate runs the toolchain's `leanchecker` (Lean v4.28 and later), which reads the compiled `.olean` of each module and adds every constant to the environment through the kernel again, with no elaborator, tactic or `unsafe` code in between; it trusts Mathlib and StatsMLlib as compiled, and the kernel that re-checks is the one of the same toolchain (not an external verifier). Negative test the same day, on a scratch package outside the repository with the same toolchain: a module that adds `false : False` with `addDeclCore (doCheck := false)` passes `lake build` and fails the gate with `(kernel) declaration type mismatch, 'false' has type Prop but it is expected to have type False`, exit 1; `#print axioms` of such a declaration would print nothing, so gate 2 does not see it (`scripts/README.md`) |
| `python3 scripts/check_core_imports.py` (gate 8, added 2026-09-08; widened the same day, F28) | under 1 s; reads the import lines of every module, no Lean | `core: 95 modules, 46210 lines (75 ours, 42494 lines; 20 vendored)`, `EXT: 177 modules, 91687 lines (157 ours, 87971 lines; 20 vendored)` (92 modules, 41403 lines, 72 ours, and 174 EXT modules, 86880 lines, 154 ours, at `0c52d54`, before the X8 count; 81 modules, 37277 lines, 61 ours, and 163 EXT modules, 82754 lines, 143 ours, at `aa27519`/`baf31a4`, before Stage 3; 37275 and 82748 lines at `da974e9`, before the docstring reflow of 2026-09-10; 81 core modules, 37825 lines, 61 ours, and 163 EXT modules, 83298 lines, at `d584cc6`, before the F36 to F38 dedup; 82 core modules, 37878 lines, 62 ours, and 164 EXT modules, 83351 lines, at `304bffc`, before F35 retired `RMT/General/Statements.lean`; 57 core modules, 24811 lines, 37 ours, and 139 EXT modules, 70284 lines, at `9c33ae6`; 70292 at `16c100c`; the 25 new core modules are the general chain `RMT/General/` and the three `Prob/` files it uses, 2026-09-09), five `entry theorem NAME: present in FILE` lines (`singleTableLaw_of_gaussian` in `RMT/Full.lean`, `heteroLaw_of_gaussian` in `RMT/Het/Sup.lean`, `tableLawR_of_gaussian_rk` in `RankR/RMT/TableLawGaussian.lean`, `heteroLawR_of_gaussian` in `RankR/Het/Sup.lean`, `singleWeightLaw_of_gaussian` in `RankR/SingleWeight/Het/Sup.lean`), `check_core_imports: PASS`, exit 0. Two independent rules: a core module may import only a core module (unchanged since 2026-09-08 morning); an EXT module may import only another EXT module, where EXT is the core plus the four Gaussian discharge families minus their three `Sup` files plus 21 definition and lemma modules, so a core violation and an EXT violation can differ (`--list` prints the core, then the EXT modules that are not core; `--boundary` prints the imports of the three `Sup` files, then of every other non-EXT module that imports an EXT module; `scripts/README.md`) |

The 12 `NO VERBATIM MATCH` lines are documented abbreviations, not mismatches. The script's own
legend states the rule: `NO VERBATIM MATCH` is expected for a signature the document
abbreviates on purpose and is a defect only when the abbreviation is undocumented, while
`NAME NOT IN TREE` is always a defect (`scripts/README.md`). Of the 12, ten are model
structures or definitions in `THEOREMS.md` section 1 where only implicit size and space binders
are dropped, and two are the Gaussian blocks of section 6.4, abridged inside the conclusion at
elaborated proof terms; one further block, in section 6.6, writes `topEigMat …` for a repeated
argument pair and the checker accepts it. No binder is dropped from the source text of any
theorem signature in that document, and a binder-by-binder comparison of the 75 declarations
quoted at commit `dfef50f` found 62 byte-identical and 13 differing only in these documented
ways (`THEOREMS.md` 0). The script reported 71 and the hand comparison 75 there because the
script skips a quoted block whose signature is shorter than 25 characters and a name it has
already compared (`scripts/check_theorems_sigs.py:23`, `if len(s)<25 or name in seen:
continue`). The 16 blocks added on 2026-09-02 all match verbatim: 87, 75, 12. The 4 blocks
of the L2 and L3 additions (sections 3.3 and 5.3) also match verbatim: 91, 79, 12.
One binder is present in the elaborated type and absent from the source text of 21 quoted
rank-one declarations: the section-level instance `[NeZero M]` (7.20, `THEOREMS.md` section 0
item 2). A source-text comparison cannot see it; the `#check` record
`notes/archive/audit_external_statement_nezero_check_2026-09-02.txt` does.

**One command for the build and four of the seven gates.** `AUDIT_PACK_JOBS=6 scripts/make_audit_pack_header.sh`
prints the UTC start time, the host, the git HEAD and worktree state, the toolchain, a full
capped build, `check_sorries.sh`, `check_axioms.sh`, `check_root_imports.py`,
`check_theorems_sigs.py`, and the UTC end time. It always exits 0 by design, so read the
recorded exit codes: a non-zero code there means the pack must not claim that the library
builds (`scripts/README.md`).

**A five-line probe** an auditor can paste into a Lean file inside the project
(`THEOREMS.md` 8):

```lean
import StackedSVD
#print axioms StackedSVD.SpikedModel.singleTableLaw_of_gaussian
#print axioms StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian
#print axioms StackedSVD.MultiTableModel.thm_svdstack_weighted_gaussian
#print axioms StackedSVD.UnalignedModelR.thm_rank_r_stacksvd_gaussian
```

Each prints `[propext, Classical.choice, Quot.sound]` (`TECHNICAL.md`). The last one was rerun
today: `lake env lean` exit 0, axioms `[propext, Classical.choice, Quot.sound]` (coordinator,
2026-09-02).

`scripts/AuditSignatures.lean` prints the elaborated signature and the axiom list of every
headline theorem, with the ambient variables and instances that a `section` hides in the
source, and it now covers the rank-`r` endpoint and the three satisfiability witnesses. Run it
with `lake env lean scripts/AuditSignatures.lean` from `lean/`. Its witness part needs
the next root build before it resolves, so it is described here and not claimed as run
(coordinator, 2026-09-02; `scripts/README.md`).

**Toolchain.** `leanprover/lean4:v4.33.0` (`lean/lean-toolchain`), Mathlib at rev
`v4.33.0` and StatsMLlib at rev `37286c3d2c7e17642fbe26988be40f3ff982f944`
(`lean/lakefile.toml`). Weyl, Davis-Kahan, Courant-Fischer and Gaussian concentration
come from StatsMLlib, not from this development (`CLAUDE.md`).

**The six numeric scripts are not gates.** Each recomputes a closed form outside Lean and
asserts it, and each prints its seed (`scripts/README.md`):
`check_port_stacksvd_subspace.py` (seconds to about 3 minutes, `OVERALL: PASS`),
`check_rayleigh_bound.py` (seconds, `RESULT: PASS`), `check_remark_examples.py` (seconds,
`RESULT: PASS`), `check_edge_count.py` (368 s), `check_singleweight.py` (seed 20260905, about
2 minutes on 2 cores, prints its rows and no PASS bar) and `check_stacksvd_orthogonality.py`
(seed 20260905, seconds, no PASS bar).

**Honest note on `check_edge_count.py`: it exits 1, by the design of its bar.** The bar at
`scripts/numeric/check_edge_count.py:225` is an exact empirical rate of `1.000`, which no finite `d`
reaches near the edge; a statement that holds with probability going to 1 gives 195 of 200
draws, that is `0.975`, and fails. All 32 failing rows are C2 or C2w rows near the edge, rows
C0, C1, C3, C4, C5 pass everywhere, the failing rates rise with `d` (for example `0.975`,
`0.990`, `1.000` at `d = 400, 800, 1600`), and no row contradicts an asymptotic claim. No Lean
statement changes under any repair option, because C2 is a numeric witness for a plan, not a
statement the tree quotes. Full analysis with the 32 rows:
`notes/archive/audit_numeric_edge_count_2026-09-02.md`; full output:
`notes/archive/check_edge_count_2026-09-02.log`.

---

## 10. Evidence trail

Independent audits, newest first. "Must-fix" is the count the note itself reports.

| Date | What it covered | Verdict | Must-fix | Note |
|---|---|---|---|---|
| 2026-09-03 | External audit of the rank-one Gaussian result on the packet `notes/audit_packets/rank1_2026-09-03/` (commit `dd64a3c`), independent, a proof-route-first protocol | 2 fatal (the paper's existential statements, `prop:binarystacksvd_inadmissable` and the two remarks, had no Lean theorem with the model built inside; closed by `Existence.lean`), 7 material (scope of the claim, the commit mismatch of the packet, the false `_inner` sentence, three disclosed-scope items), 4 minor (the file-layering sentence of section 4, false for the Het chain; 7.16 attached to two wrong results; a duplicated packet row; the empty transcript appendix). The audit's own proof route agrees with the implemented one at every step and found no shorter verified route | 2, closed | `notes/archive/audit_external_rank1_2026-09-03.md`, response `notes/archive/audit_external_rank1_response_2026-09-03.md` |
| 2026-09-03 | L5 (the Gaussian marginalization behind `app:wstacksvd_mle`), coordinator's adversarial audit of the seven new files, with a numeric check | pass: the random-effects law matches `SpikedModel.X`, the explicit square root and the density constant verified numerically (seed 20260903, `A Aᵀ = Σ` to `1e-16`, grid integral of the density `0.9999995`, clause 2 to `4e-15`), the theorem is not vacuous | 0 | `notes/archive/L5_mle_marginal.md` (section "Adversarial audit"), `notes/archive/audit_L5_numeric_2026-09-03.py` |
| 2026-09-02 | Cleanup pass 2 (commit `c858dba`), adversarial | loses no mathematics; 24 deletions all redirect to a survivor, 8 hypothesis drops sound | 3 (record and docstrings), 7 should-fix | `notes/archive/audit_cleanup_pass2_2026-09-02.md` |
| 2026-09-02 | Cleanup pass 1 (commit `9d727cc`), adversarial | loses no mathematics; 14 deletions, 5 dropped hypotheses sound | 3 (record and docstrings), 7 should-fix | `notes/archive/audit_cleanup_pass1_2026-09-02.md` |
| 2026-09-02 | `TECHNICAL.md` and `THEOREMS.md` as auditor documents, ten checks | ready for an external auditor; 74 edits in five files; nothing found was a false statement about the mathematics | 0 | `notes/archive/audit_docs_2026-09-02.md` |
| 2026-09-02 | `check_edge_count.py` exit 1 | not a failed claim; the bar is exact 1.000 | 0 | `notes/archive/audit_numeric_edge_count_2026-09-02.md` |
| 2026-09-02 | Track E wave 2b (E7, E8bd, E8a, E8b, E9), independent | all five files pass; instantiated the final theorem on a concrete Gaussian model and matched `gammaR` to `theory_pred.R` to `1.0e-13` | 0, 9 should-fix (one a scope note) | `notes/archive/audit_independent_TrackE_wave2b_2026-09-02.md` |
| 2026-09-02 | Track E wave 2a (E3, E4, E6), independent | all three files pass | 0, 5 should-fix | `notes/archive/audit_independent_TrackE_wave2a_2026-09-02.md` |
| 2026-09-02 | Track E wave 1 (E0, E1, E2, E5, E10), independent | all six files pass | 0, 6 should-fix | `notes/archive/audit_independent_TrackE_wave1_2026-09-02.md` |
| 2026-09-02 | Track E plan, before any proof | accept with corrections; the target statement was wrong as written and was restated | see note | `notes/archive/audit_rankr_plan_E_2026-09-02.md` |
| 2026-09-02 | Track D (the rank-`r` aligned appendix) | pass with findings; 0 fatal | 0, 9 minor | `notes/archive/audit_independent_trackD_2026-09-02.md` |
| 2026-09-02 | Track C (Gaussian `TableLawR` at general `r_i`) | pass with findings; 0 fatal | 0, 7 minor | `notes/archive/audit_independent_trackC_2026-09-02.md` |
| 2026-09-02 | Track A (rank-`r` stacksvd for Gaussian noise) | pass with findings; 0 fatal | 1 (scope documentation, since addressed), 4 minor | `notes/archive/audit_independent_trackA_2026-09-02.md` |
| 2026-09-02 | Track C plan, before any proof | accept with corrections; one route dropped from 900 to 1400 lines to 200 to 350 | see note | `notes/archive/audit_rankr_plan_C_2026-09-02.md` |
| 2026-09-02 | Stage D33 rank-1 batch (P1, P2, P3, MLE, remarks), independent | pass with findings; every theorem with more than two hypotheses instantiated on a concrete Gaussian model, so no hypothesis set is contradictory | 0, 5 minor | `notes/archive/audit_independent_D33_2026-09-02.md` |
| 2026-09-01 | Track B (general-`r_i` svdstack), independent | on track; 0 blockers | 0, 3 weakenings, 8 cosmetics | `notes/archive/audit_independent_trackB_2026-09-02.md` |
| 2026-09-01 | Third external review of audit pack v3, plus the verdicts | five holds, two overstated, several superseded | see note | `notes/archive/external_audit_v3_review_2026-09-01.md`, `..._verdicts_2026-09-01.md` |
| 2026-09-01 | Track A and B plans, before any proof | both plans survive with corrections | see notes | `notes/archive/audit_rankr_plan_A_2026-09-01.md`, `..._B_2026-09-01.md` |
| 2026-09-01 | `MLE.lean` and `Remarks.lean` against the paper | no vacuity, no mismatch, no hidden contradiction | 0, 4 low severity | `notes/archive/audit_mle_remarks_2026-09-01.md` |
| 2026-09-01 | The Gaussian discharge of `HeteroLaw`, items H0 to H11 (11 files, 11989 lines) | see note | see note | `notes/archive/audit_het_chain_2026-09-01.md` |
| 2026-08-31 | Second external review of pack v2, plus the verdicts | main new claim credible at the level of the quoted statements | see note | `notes/archive/external_audit_v2_review_2026-08-31.md`, `..._verdicts_2026-08-31.md` |
| 2026-08-31 | Whole project against the paper and `notes/archive/PLAN_2026-08-29.md` | see note (findings 3, 4 drove decisions D20, D21) | see note | `notes/archive/audit_global_2026-08-31.md` |
| 2026-08-31 | Mechanical failure modes (unused hypotheses, junk values) | 3 of 8 headline theorems carry an unused hypothesis, each documented | see note | `notes/archive/audit_mechanical_2026-08-31.md` |
| 2026-08-31 | The five pieces that no audit had covered | faithful, with the recorded additions | see note | `notes/archive/audit_unaudited_2026-08-31.md` |
| 2026-08-30 | First external audit (statement level, proof bodies removed), plus the verdicts | the headline finding was correct and changed the plan to the sharp Sudakov-Fernique route | see notes | `notes/archive/external_audit_2026-08-30.md`, `..._verdicts_2026-08-30.md` |
| 2026-08-30 | The Gaussian proof chain of `prop:single_table` | the chain is sound and not vacuous | see note | `notes/archive/audit_l2_chain_2026-08-30.md` |
| 2026-08-29 | Every closed form, recomputed against the Python and R reference code (seed 20260829) | worst difference `3.3e-16`; one disagreement, the degenerate `A_β = I` case, which Lean excludes | see note | `notes/archive/audit_numeric_2026-08-29.md` |
| 2026-08-29 | Scope and framing of the first Lean files | 24 faithful, 7 stronger hypothesis, 2 weaker conclusion, 1 different object, 5 not in paper; two findings blocked Layer 1 and were fixed | see note | `notes/archive/audit_scope_2026-08-29.md` |
| 2026-08-29 | The external audit of the plan itself, which drove the v2 interface | directionally sound, not the best Lean architecture; four changes adopted | see note | `notes/archive/lean_audit.txt` |

`notes/archive/` also holds the four audit packs (`AUDIT_PACK.md`, `AUDIT_PACK_V2.md`,
`AUDIT_PACK_V3.md`, `AUDIT_PACK_HETEROEDGE.md`) sent to earlier external reviewers, and
`notes/archive/AUDIT_PACK_V4.md` is the current one (`CLAUDE.md`, repository listing).

---

## 11. Decisions and open questions

**Decisions taken without the user.** `notes/FLAGGED.md` is the register. Its table has 26
numbered rows, `D1` to `D24` plus `D26` and `D27`; there is no `D25` (grep of the file, this
session). Five rows, `D19`, `D23`, `D24`, `D26` and `D27`, record the user's own answers rather
than a decision taken alone. Seven further blocks, `D28` to `D34`, carry decisions in prose,
and the `D32` block carries a numbered sub-list of 17 items, several with their own
sub-decisions.

The nine decisions once marked "revisit: yes" (D5, D6, D12, D13, D17, D18, D20, D21, D22)
were all closed by the user on 2026-09-05; no decision waits for the user (`notes/FLAGGED.md`,
`notes/SERVER_TODO.md`).

The ones a reader should look at first, because each one moves a statement:

| # | Decision | Why it matters |
|---|---|---|
| D6 | `delocUniform` proved from right-rotation invariance of the Gaussian model; `ResolventLimits` loses six `w`-fields | changed the shape of the Layer 2 interface |
| D12, D22 | weighted stackSVD proceeds over `HeteroLaw`, then route A (`RMT/Het/`) discharges it | the hypothesis is now proved, so the decision is closed by `heteroLaw_of_gaussian` |
| D13 | rank-`r` performance is `tr(Bᵀ specInvTop A r B)` | this is the trace functional of scope item 7.8 |
| D14, D16 | `Rank(∑ R_i R_iᵀ) = r` dropped from the model; `rank B_R = r` chosen as the replacement hypothesis after the weaker form was found false by counterexample | changes the hypothesis of `thm:gen_rank_weight_svdstak` |
| D17, D18 | `SubspaceLaw` keeps `align` only; `topGap` dropped as unsatisfiable when `∑_i n_i N < r ≤ d N` | changes a structure the reader will check against the paper |
| D21 | `HeteroLaw.align` stated at `L(w)`, one step past the cited display | scope item 7.13 |
| D-A1 (item 10) | document the ordered-class narrowing now, move to an order-free model later if the user wants | scope items 7.2, 7.3 |
| D23 | two appendix results go to the end of the queue | scope item 7.5 |

**Open questions.** `Q1` to `Q16` (`notes/FLAGGED.md`, `THEOREMS.md` 7.3), one line each:

| # | Question | State |
|---|---|---|
| Q1 | The R package `theory_pred.R` fails at `M = 1` above threshold and at `∑θ⁴/c = 1` | outside the repository, not touched |
| Q2 | `A_β = I`: the Python and R code disagree and the paper leaves it open | Lean excludes it through `hthr` |
| Q3 | Layer 1 is Gaussian-only through `JointGaussianNoise`; `assum:general_noise` is wider | keep, user decision 2026-08-29; since 2026-09-02 the rank-one SVDstack family takes `IndepNoise` (L4); since 2026-09-08 `thm_theta_est` takes `ThetaEstLaw` (F31), so no rank-one Layer 1 theorem keeps `JointGaussianNoise`; since 2026-09-08 (Stage 0 of `notes/NONGAUSSIAN_SCOPE.md`) `ThetaEstLaw` is also discharged for any fixed law with mean 0, variance 1 and finite fourth moment (`thetaEstLaw_of_general`, `thm_theta_est_general`), the first non-Gaussian discharge |
| Q4 | Prove a block-variance Marchenko-Pastur law, or stop at a second hypothesis structure | closed: `heteroLaw_of_gaussian` |
| Q5 | The `M = 3` remark example sits on the threshold | now `notes/paper_edits.md` E6 |
| Q6 | `thm:simple_thm1` is cor. 1, not cor. 2 | fixed in `notes/README.md` |
| Q7 | Five rank-`r` design decisions | answered; item (5) closed by the Gaussian rank-`r` law |
| Q8 | The removal remark at line 901 is wrong in general | not formalized, scope item 7.5 |
| Q9 | The paper writes `prop:stacksvd_general` where it means `prop:stacksvd_subspace` (line 842) | label typo, for the authors |
| Q10 | Four appendix results had no Lean | two are now proved (`thm:rank_r_stacksvd`, `thm:rank_r_svdstack`); the other two now have Lean at Layer 1 (`prop:gen_rank_stacksvd_singleweight`, `prop:singleweight_suboptimality`, 2026-09-05), with the Gaussian discharge of `SingleWeightLaw` open |
| Q11, Q12 | The exact heteroscedastic edge | closed by `heteroEdge_of_gaussian` |
| Q13 | External audit verdicts and the sharp edge route | closed by D28 |
| Q14 | The second external review of the audit pack | verdicts recorded |
| Q15 | The MLE and remark statements | answered "OK both", 2026-09-01 |
| Q16 item 1 | Keep the ordered rank-`r` class, or add a permutation `R_i` route for the paper's worked example | **open**; Claude recommends keeping the class |
| Q16 item 2 | The component-clause statement form of `thm:rank_r_svdstack` (canonical frame, `hS : StrictAnti`) | **closed**: user OK 2026-09-05 ("prove simplest rank r (S_I sorted)"), `notes/FLAGGED.md` Q16 |
| Q16 last | Whether `notes/paper_edits.md` E1 to E11 go to the paper before or after this audit | **closed** 2026-09-05: after; the audit reads the paper as it is, with E1 to E11 as the list of clarifications the authors will make |

---

## 12. Repository map

From `TECHNICAL.md` (layout); since 2026-09-10 `CLAUDE.md` holds only the working rules and a
ten-line map of the top level, and every folder has a README.

**Top-level files** (nine entries since the tidiness pass of 2026-09-10, commits `baf31a4` and
`66be99d`). `README.md` (the plain-language entry point: a path per reader, the result table,
the caveats, how to check it, the folder map), `GETTING_STARTED.md` (the
coauthor's install-and-check guide, at the top level since 2026-09-10), `CLAUDE.md` (the
working rules, 104 lines), `CITATION.cff` (since 2026-09-10), `LICENSE`, and the four folders
`docs/`, `lean/`, `notes/`, `scripts/`. There is no continuous integration (a user decision
of 2026-09-10; the gates run by hand on the server). **`docs/`** holds `TECHNICAL.md` (the
detailed paper-to-theorem table, the notation table, the full file tree and the verification
section; the top-level README until 2026-09-03), `THEOREMS.md` (the auditor annex),
`notes/paper_edits.md` (findings E1 to E12), `notes/NOTE_FOR_COAUTHORS.md` (one page), `SORRIES.md` (the
`sorry` ledger, empty since 2026-09-10 evening, one tracked row for part of that day (Stage
3's Furedi-Komlos count), empty before that; read by
`scripts/check_sorries.sh`; at the top level until 2026-09-10),
`paper/main_paper.tex` (a read-only snapshot of the paper, taken 2026-08-29; at the top level
until 2026-09-10), a `README.md` that names each, and this file. The plan of
record for 2026-08-29 to 2026-08-30 is `notes/archive/PLAN_2026-08-29.md` (kept for the decision
trail, not the current state).

**`lean/`** is the lake project (`StackedSVD/` until the rename of 2026-09-10). `lean/StackedSVD/` is the source. Measured on the
tree of the build of record: **108508 lines in 210 files ours, plus 3716 vendored lines in 20
files** (`wc -l` over the tree, coordinator, 2026-09-10 16:03 EDT, on the tree of
commit `416ab9d`; 103701 lines in 207 files at `0c52d54` (the Stage 3 commit before the X8
count, 2026-09-10); 99510 lines in 196 files at the rename commit
`baf31a4`, the same count as at 07:30 EDT after the docstring reflow of the path fix; 99504 lines in 196 files at 2026-09-09 23:29 EDT, after the F36 to F38 dedup
of `RMT/General/`; 100054 lines in 196 files at 21:50 EDT, after F35 retired
`RMT/General/Statements.lean`; 100105 lines in 197 files at 20:35 EDT, after the close of the
non-Gaussian Stage 1; 99464 lines in 195 files at 20:00 EDT, after wave 5 of the same stage;
86852 lines in 171 files at 2026-09-08 23:40 EDT, after F32; 86894 lines in
171 files at 21:30 EDT, after F28; 86623 lines in 165 files at 16:05 EDT, after Stage 0; 85630 lines
in 162 files at 2026-09-07 21:50 EDT, after `Main.lean` and F30;
84840 lines in 161 files at 19:10 EDT, after the linter campaign;
82395 lines in 154 files at `9c82538`, 2026-09-05 17:00 EDT, after Track G; 78007 lines
in 146 files at `947b453`; 77552 lines in 145
files at `482309c`; 76561 lines in 139 files at `697da9c`; 76134 lines in 137
files at `a5decd4`; 76011 lines in 137
files at `04fd96b`; 75744 lines
in 136 files at `e2c76c5`, 74563 lines
in 129 files at `ea88ba3`, 74537 at `e3a2912`, 74498 at `ea1763d`, 74268 lines in 127 files
at `e5cff73`).

| Directory | Files | Content |
|---|---|---|
| root of the source | 16 | the rank-one model (`Defs.lean`), `Spectral.lean`, the hypothesis structures `RMT.lean` and `StackSVDWeighted.lean` (its Layer 1 consumer moved to `StackSVD/Weighted.lean`, F28), `StackSVD.lean` (likewise, to `StackSVD/Main.lean`), `Scalars.lean`, `Secular.lean`, `ThetaEst.lean`, `MLE.lean`, `MLEConverse.lean`, `Remarks.lean`, `RemarksUniform.lean`, `StrictFacades.lean`, `SVDStack.lean`, `Sat.lean`, `Existence.lean` (section 6) |
| `StackSVD/` | 2 | new 2026-09-08 (F28): `Main.lean` (`prop:stacksvd_general`, `thm:simple_thm1` stacksvd half) and `Weighted.lean` (`cor.2` binary, `thm:stacksvd_weighted`, `thm:stacksvd_binary_optimal_svd_stack`), the Layer 1 theorems moved out of `StackSVD.lean` and `StackSVDWeighted.lean` |
| `SVDStack/` | 11 | `A_β`, the svdstack estimator, `thm:svd_stack_general`, `thm:svdstack_weighted`, the uniform Rayleigh bound, `thm:simple_thm1`, the signed form of `lem:entrywise_conv_eigenvec` (`EntrywiseSigned.lean`), the svdstack clauses of `prop:binarystacksvd_inadmissable` (`Inad.lean`), the finite-`d` form of `lem:delocalization` (`DelocFinite.lean`, F1); since F28, `DelocDir.lean`, the consumer-free half of `Gram.lean` |
| `LinAlg/` | 10 | Weyl, Davis-Kahan, top-`r` and index projector perturbation, Ky Fan, frames |
| `Prob/` | 6 | convergence in probability, polynomial null sets, Gaussian matrix laws and adapters; since 2026-09-03 `WithDensityPi.lean` (densities of product measures and of linear images) and `GaussianDensity.lean` (the density of `N(0, S)` on `Fin d → ℝ`) |
| `MLEMarginal/` | 5 | the random-effects model of `app:wstacksvd_mle` and its density: `Defs`, `RowLaw`, `TableLaw`, `LogLik`, `Main` (`thm_wstacksvd_mle_marginal`) |
| `RMT/` | 19 | Layer 2: the Marchenko-Pastur chain, `singleTableLaw_of_gaussian` |
| `RMT/Het/` | 14 | Layer 2 heteroscedastic: `heteroLaw_of_gaussian`, `heteroEdge_of_gaussian`, `thm_stacksvd_weighted_gaussian`, `prop_dominance_gaussian` |
| `RMT/General/Edge/` | 14 | Stage 3 (2026-09-10): the sharp upper edge of the noise Gram matrix at four moments by the moment method with a truncation (route C of `notes/stage3_edge.md`), `opNorm_sq_edge_of_general`, `lamMax_W0_edge_of_general`, `singleTableLaw_of_moments`; the Furedi-Komlos shape count `cellCard_mul_le_pos` (`CountBound.lean`), proved 2026-09-10 evening by the X8 plan (`notes/x8_plan.md`) |
| `RankR/` | 25 | Section 7 statements and Layer 1: the two models, `TableLawR`, `SubspaceLaw`, `HeteroLawR`, the four Section 7 results, asymptotic orthogonality of the svdstack columns (`AlignedOrth.lean`, F2); since F28, `Unweighted.lean`, `WeightedMain.lean`, `SubspaceMain.lean` hold the Layer 1 theorems moved out of `Defs.lean`, `Weighted.lean` and `Subspace.lean` |
| `RankR/RMT/` | 25 | Gaussian `TableLawR` at general `r_i` and Gaussian `SubspaceLawG` |
| `RankR/Het/` | 13 | Track E: `heteroLawR_of_gaussian`, `thm_rank_r_stacksvd_gaussian`; `Sub` (F8, the sub-model that drops the zero-weight tables) |
| `RankR/SingleWeight/` | 9 | one weight per table (Appendix D): `SingleWeightLaw`, `prop_gen_rank_stacksvd_singleweight`, `prop_singleweight_suboptimality_of_law`, `prop_singleweight_suboptimality_gaussian` (section 6.8, 6.9; `Tie.lean`, `Regimes.lean`, `Suboptimality.lean` are F18b) |
| `RankR/SingleWeight/Het/` | 8 | Track G: `singleWeightLaw_of_gaussian` and its facades |
| `Vendor/COLT83/` | 20 | Stein, Gaussian interpolation, Sudakov-Fernique, Borell-TIS (backport of `2d84602`) |

`lean/StackedSVD.lean` imports every module.

**`scripts/`** holds the eight gates (`check_sorries.sh`, `check_axioms.sh` with its driver
`AxiomAudit.lean`, the coverage gate `check_root_imports.py`, the document gate
`check_theorems_sigs.py`, the layering gate `check_layering.sh` with its driver
`LayerAudit.lean`, section 4, and the paper-edit gate `check_paper_edits.py`, which applies the
exact replacements of `notes/paper_edits.md` to a copy of the paper snapshot in memory and checks
that every "Paper, as amended" block of `THEOREMS.md` is verbatim in the result, the kernel gate `check_kernel.sh`, which replays every module through the kernel with the toolchain's `leanchecker`, and the core gate `check_core_imports.py`, which checks that the 92 modules of the Gaussian RMT core (57 at its creation on 2026-09-08; 81 before Stage 3, 2026-09-10) import no paper-specific module), the document path check `check_paths.py` (2026-09-10, not a gate of record: every path the reader-facing documents name exists), the one-command header `make_audit_pack_header.sh` that runs
the build and the first four in one pass, the signature dump `AuditSignatures.lean`, the
six numeric checks, the packet builder `make_audit_packet_rank1.py`, the exploration
scripts of `scripts/numeric/` (the top-level `scratchpad/` until 2026-09-09), and
`scripts/README.md` documenting all of them.

**`notes/`** is the lab notebook of the development repository; it is not part of this
public copy, and the `notes/` paths that this document and the Lean docstrings cite name
its files (`scripts/release_withheld.txt` lists them). Living files: `PROGRESS.md` (the log; its "Where we are"
block is the fastest catch-up), `SERVER_TODO.md` (per-module status and what comes next),
`FLAGGED.md` (decisions and questions), `INTERFACES.md` (the fixed names every note and agent
uses), `README.md` (every paper label with its formalization class A, B or C; `00_inventory.md`
until 2026-09-10, and the translation of the paths that archived material still uses),
`FOLLOWUP_LIST.md` and `REVISION_LIST.md` (the follow-up items and the paper revisions that
wait), `RANK_R_PLAN.md`, `NONGAUSSIAN_SCOPE.md` and `stage3_edge.md` (the non-Gaussian
extension and its deferred edge stage), `REPO_CLEANUP_PLAN.md` (2026-09-10),
`BUILD_HISTORY.md` (every earlier build of record, moved out of `CLAUDE.md` and
`TECHNICAL.md` on 2026-09-10), `SERVER_SETUP.md` (the server and laptop setup, out of
`CLAUDE.md` since 2026-09-10), `discoveries/` and `audit_packets/` (each with a README).
`notes/archive/` holds completed material: per-theorem review notes, earlier audits and
audit packs, external reviews, campaign plans, and `agent_reports/` (one report per proof
agent, 197 files; the 93 files of `notes/agent_reports/` folded in on 2026-09-10), listed in
its `INDEX.md`.


---

## Appendix A. The build and gate transcript

Verbatim output of one command, `AUDIT_PACK_JOBS=6 scripts/make_audit_pack_header.sh`, run
2026-09-10 16:18:21 to 16:19:22 EDT (20:18:21 to 20:19:22 UTC, 61 s: 0 jobs rebuilt, 5
replayed, every other job up to date from the 15:58:20 to 16:01:47 EDT root build) on
`416ab9d` (`Stage 3 complete: the Furedi-Komlos count proved (X8 plan), 0 sorry`), with a
clean worktree (0 files modified or untracked at header time), so the transcript audits
exactly the committed tree. It prints the git
state, the toolchain, the full build and four of the eight gates; the fifth,
`scripts/check_layering.sh`, the sixth, `scripts/check_paper_edits.py`, the seventh,
`scripts/check_kernel.sh`, and the eighth, `scripts/check_core_imports.py`, are run
separately (section 9, their rows). The script always exits 0 by design, so read the
recorded exit codes, not its own. The two numbers after `declarations audited` are the
auditor's `<reported> <skipped>` counts (declarations reported, compiler-generated
auxiliaries skipped; `scripts/AxiomAudit.lean`). The `check_axioms` WARNING line is the script's note that
the header run itself did not compile (`CHECK_AXIOMS_BUILD=0`): the gate reads the oleans of
the 15:58:20 to 16:01:47 EDT build, on the committed `416ab9d` tree; its output shows 0
`sorryAx` and no tracked sorry. The transcripts of `aa27519` (2026-09-10 13:04 EDT),
`66be99d` (09:33 EDT), `6633055` (08:23 EDT),
`b8771cb` (2026-09-09 23:41 EDT), `d584cc6` (2026-09-09 21:51 EDT), `304bffc` (2026-09-09
20:43 EDT), `9c33ae6` (2026-09-08 23:35 EDT), `16c100c` (21:58 EDT), `9fdcbfa` (16:50 EDT),
`8736b2a` (08:32 EDT), `290d485` (2026-09-07 22:32 EDT), `a32907a` (2026-09-07 21:06 EDT),
`9c82538` (2026-09-05 19:36 EDT), `a2cbbb9` (16:58 EDT), `947b453` (13:52 EDT), `697da9c`
(10:41 EDT) and `a5decd4` (2026-09-03) that stood here before are in the git history of
this file.

```
generated-by : scripts/make_audit_pack_header.sh
start (UTC)  : 2026-09-10 20:18:21 UTC
host         : (the development server)

--- git ---
HEAD         : 416ab9dfa3196b4dd7a04de77024ad6e7312e7ba
log -1       : 416ab9d Stage 3 complete: the Furedi-Komlos count proved (X8 plan), 0 sorry
worktree     : 0 file(s) modified or untracked at header time

--- toolchain ---
lean-toolchain: leanprover/lean4:v4.33.0
lake version  : Lake version 5.0.0-src+d8b1897 (Lean version 4.33.0)

--- build: lake-build-capped.sh 6 ---
last line    : Build completed successfully (8978 jobs).
exit code    : 0
log lines    : 144
errors       : 0
jobs rebuilt : 0
jobs replayed: 5  (cached traces re-checked, not recompiled)

--- gate: scripts/check_sorries.sh ---
check_sorries: OK - 0 sorry sites in 0 declarations, 0 rows, all matched.
exit code    : 0

--- gate: scripts/check_axioms.sh ---
check_axioms: WARNING - no build. The audit reads the oleans in .lake/build.
  newest source: 2026-09-10 15:58:13  <repository>/lean/StackedSVD/RMT/General/Edge/CountBound.lean
  oldest olean : 2026-09-02 13:47:43  <repository>/lean/.lake/build/lib/lean/StackedSVD/Vendor/COLT83/Mathlib/Probability/SteinIdentity.olean
check_axioms: 5631 declarations audited (5631 3297).
check_axioms: OK - no axiom outside {propext, Classical.choice, Quot.sound} and no untracked sorry.
exit code    : 0

--- gate: scripts/check_root_imports.py ---
check_root_imports: 230 modules (210 ours, 20 vendored); 229 reachable from lean/StackedSVD.lean.
  vendored and unreachable (not gated, see the module docstring):
    StackedSVD.Vendor.COLT83.Axioms
check_root_imports: OK - every module of the package is audited.
exit code    : 0

--- gate: scripts/check_theorems_sigs.py docs/THEOREMS.md ---
NO VERBATIM MATCH: SpikedModelR  in: ['RankR/General.lean']
NO VERBATIM MATCH: UnalignedModelR  in: ['RankR/GramR.lean', 'RankR/General.lean']
NO VERBATIM MATCH: TendstoInProb  in: ['ThetaEst.lean', 'Defs.lean']
NO VERBATIM MATCH: overlap  in: ['Defs.lean']
NO VERBATIM MATCH: Abeta  in: ['SVDStack/Defs.lean']
NO VERBATIM MATCH: stackSVDLimit  in: ['StackSVD.lean']
NO VERBATIM MATCH: svdstackLimit  in: ['SVDStack/Defs.lean']
NO VERBATIM MATCH: Sval  in: ['SVDStack/Defs.lean']
NO VERBATIM MATCH: optW  in: ['SVDStack/Defs.lean']
NO VERBATIM MATCH: svdstackLimitOpt  in: ['SVDStack/Defs.lean']
NO VERBATIM MATCH: thm_gen_rank_weight_svdstak_general_r_paper_gaussian  in: ['RankR/GeneralGaussian.lean']
NO VERBATIM MATCH: thm_gen_rank_weight_svdstak_general_r_full_gaussian  in: ['RankR/GeneralGaussian.lean']
signatures 148, verbatim 136, no verbatim match 12
A 'NO VERBATIM MATCH' line means the quoted signature does not appear character for
character in the tree. That is expected for a signature the document abbreviates on
purpose; docs/THEOREMS.md says which ones and why. It is a defect only when the name is
there and the abbreviation is not documented. 'NAME NOT IN TREE' is always a defect.
exit code    : 0

end (UTC)    : 2026-09-10 20:19:22 UTC
```
