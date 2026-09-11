# scripts/ - mechanical gates before a commit

The working rules of the development repository (`CLAUDE.md`, not part of this copy)
require the gates to pass before a commit; the first two are the
minimum, the eight together are the build of record.

```sh
scripts/check_sorries.sh && scripts/check_axioms.sh
```

Run them from anywhere. Each script finds the repository root from its own path.

## check_sorries.sh

It compares the `docs/SORRIES.md` table with the `sorry` sites in the tree.

1. It reads every `lean/StackedSVD/**/*.lean` file, without `Vendor/`.
2. It blanks line comments, nested block comments, doc comments, string literals and
   character literals, and keeps the newlines, so the line numbers do not move. The
   word `sorry` inside a comment or inside a message is therefore not a site.
3. It matches `sorry` as a whole token. `sorryAx`, `notsorry` and `Foo.sorry` do not
   match.
4. For each site it takes the nearest preceding `theorem`, `lemma`, `def`, `instance`,
   `abbrev`, `example`, `structure`, `class` or `inductive` line as the declaration.
5. It reads the `File` and `Declaration` columns of `docs/SORRIES.md`. The `File` cell
   `StackedSVD/<path>.lean` means `lean/StackedSVD/<path>.lean`.
6. It compares the pairs (file, last dotted component of the declaration).

Options:

| Option | Effect |
|---|---|
| (none) | gate |
| `--sites` | print `file:line<TAB>declaration` for every site |
| `--tracked` | print the declaration of every `docs/SORRIES.md` row, one per line |
| `--help` | usage |

Exit codes: `0` match or listing, `1` mismatch (the diff names the rows without a site
and the sites without a row), `2` usage or environment error.

Limits. A `sorry` inside an anonymous `instance` gets the label
`<anonymous instance line N>`, which no `docs/SORRIES.md` row can match; name the instance.
A declaration header that spans two lines (the keyword on one line, the name on the
next) is not parsed. Two `sorry` sites in one declaration count as one row.

## check_axioms.sh and AxiomAudit.lean

`AxiomAudit.lean` imports `StackedSVD`, walks every module of the `StackedSVD`
namespace, and calls `Lean.collectAxioms` on each declaration (the same function that
`#print axioms` and Mathlib's `assert_no_sorry` use). It prints one record per line;
`check_axioms.sh` parses the records and gives the verdict.

The rule: a declaration passes when its axioms are inside
`{propext, Classical.choice, Quot.sound}`. `sorryAx` is the one conditional case.

1. The declaration has a `docs/SORRIES.md` row: pass, listed as `tracked sorry`.
2. The declaration only uses another declaration that has a row: pass, listed as
   `inherits a tracked sorry`. The driver reports the origin of each `sorry` on a
   `SORRYDEPS` line, so an inherited `sorry` needs no row of its own.
3. Anything else: fail.

Any axiom outside the three also fails. `CLAUDE.md` hard rule 2 forbids `axiom`, so this
gate is how a new `axiom` gets caught.

Options and variables:

| | Effect |
|---|---|
| `--raw` | print the driver output, no verdict |
| `--verbose` | gate, and print every declaration with its axioms |
| `--help` | usage |
| `CHECK_AXIOMS_BUILD=1` | build first: `lake-build-capped.sh 12` on the server, `lake build` elsewhere (see the caveat) |
| `AXIOM_AUDIT_VENDOR=1` | audit `StackedSVD/Vendor/` too (skipped by default) |

Exit codes: `0` pass, `1` axiom violation, `2` usage, environment, build or driver error.

### CAVEAT: stale oleans

The driver reads the compiled `.lake/build` oleans, not the sources. Without
`CHECK_AXIOMS_BUILD=1` the script does not build. It prints the newest source mtime, the
oldest olean mtime, a `STALE` line when a source is newer, and a `NO OLEAN` line for
every source file that the build does not carry. Read those lines before you trust the
verdict:

- a declaration edited since the last build is audited in its old form;
- a declaration added since the last build is not audited at all;
- a declaration deleted since the last build is still audited.

`CHECK_AXIOMS_BUILD=1 scripts/check_axioms.sh` is the real gate. It refuses to start when
another `lake` process runs, because two builds on one `.lake/build` corrupt each other.
On the development server it builds with `$LEAN_TOOLS/lake-build-capped.sh 12` (never run
`lake build` uncapped there); on any other machine it runs `lake build`.

Where the scripts run: `check_axioms.sh`, `check_layering.sh`, `check_kernel.sh` and
`make_audit_pack_header.sh` read `LEAN_TOOLS`, the folder of the helper scripts of the
development server (`env.sh`, `lake-build-capped.sh`, `lean-local.sh`, the nprocs shim;
the default of this public copy is empty). They load the toolchain from `$LEAN_TOOLS/env.sh` when that file
exists and otherwise use the elan-managed `lake` on `PATH`; the timestamps come from
`python3`, so the scripts do not need GNU `find` or `date`. Tested on the server (Linux);
not yet run on macOS.

A `lake build` that runs at the same time rewrites the oleans under the driver and can
make it exit non-zero with no message. The script retries once in that case.

### Limits of the driver

It skips, and counts as skipped:

- compiler generated names: `_private` internals, `.proof_N`, `.match_N`, `.eq_N`,
  `_auxLemma.N`, every `_`-prefixed part, and names whose last component is `rec`,
  `recOn`, `casesOn`, `brecOn`, `injEq`, `noConfusion`, `mk`, `inj` and similar;
- `Decidable`, `DecidableEq` and `DecidablePred` instances.

Private declarations are audited. They appear under their user facing name with the
suffix ` (private)`.

The environment is what `import StackedSVD` pulls in, so a file that
`lean/StackedSVD.lean` does not import is never audited, even when its olean
exists. Add every new module to that root file.

## check_layering.sh and LayerAudit.lean

`LayerAudit.lean` imports six `StackedSVD` modules directly (not the root; see the file's
doc comment) and walks the transitive closure of project constants that each of the three
Layer 2 Gaussian discharges (`singleTableLaw_of_gaussian`, `heteroLaw_of_gaussian`,
`heteroEdge_of_gaussian`) uses in its statement and its proof. Two criteria flag a
constant in that closure: a fixed list of 13 Layer 1 theorem names, such as
`thm_stacksvd_weighted` and `thm_theta_est` (`FORBIDDEN` lines), and a list-free one, any
constant whose type takes a `SingleTableLaw` or a `HeteroLaw` as a hypothesis, also under a
`∀ i` (`CONSUMER` lines; the structures' own projections such as `HeteroLaw.align` are
exempt). Before the audit the driver self-tests the second criterion on five known
consumers and on the three endpoints (`SELFTEST` line; a failure is a Lean error). Both
criteria are evidence that a discharge does not depend on the proof of the theorem it
exists to feed, a question the file-level import graph cannot answer on its own
(`RMT/Het/MPhet.lean` and similar import `StackSVDWeighted.lean` only for the type
`HeteroLaw`). A negative control of 2026-09-03 (a copy of the driver with
`thm_stacksvd_weighted_gaussian` added as a fourth endpoint) gave 2 `FORBIDDEN` and 1
`CONSUMER` lines and `LAYERAUDIT_END 4 3`.

`check_layering.sh` runs the driver under a 1500 s timeout and gives the verdict: pass when
all 3 endpoints are audited and 0 forbidden hits turn up. `--raw` prints the driver output
with no verdict; `--help` usage.

Exit codes: `0` pass, `1` a forbidden constant turned up or the endpoint count is wrong, `2`
usage, environment, or Lean error. State on 2026-09-03 (oleans of the 09:02 EDT root build):
3 endpoints, closures of 1288, 1562 and 323 project constants, 97 `USES` lines over 75
distinct constants (the structure `HeteroLaw`, its constructor, and definitions and
structural lemmas of `StackSVDWeighted.lean` such as `stack`, `stackGramW`, `Scalars.Lw`,
`Scalars.secular`), 0 `FORBIDDEN`, 0 `CONSUMER`, exit 0, 30 s.

## check_root_imports.py - the coverage gate

```sh
python3 scripts/check_root_imports.py          # gate
python3 scripts/check_root_imports.py --list   # every module, with its status
```

`AxiomAudit.lean` audits what `import StackedSVD` pulls in, so a module that
`lean/StackedSVD.lean` does not reach, directly or transitively, is compiled by
`lake build` and then never audited. A `sorry` or an `axiom` in such a module passes both
gates unseen. This script computes the transitive closure from the root file and fails when a
module of ours is outside it. It reads the `public import` and `private import` forms of the
Lean 4.33 module system as well as plain `import`.

Vendored code (`StackedSVD/Vendor/`) is reported, not gated: the axiom driver skips it unless
`AXIOM_AUDIT_VENDOR=1`. One vendored module is outside the closure on purpose,
`Vendor/COLT83/Axioms.lean`, which is the vendor's own `#print axioms` driver.

State on 2026-09-02, after `StackedSVD/Sat.lean` landed: 147 modules, 127 ours, all 127 reachable.

Exit codes: `0` pass, `1` a module of ours is unreachable, `2` environment error.

## make_audit_pack_header.sh - one command for an external auditor

```sh
AUDIT_PACK_JOBS=6 scripts/make_audit_pack_header.sh
```

It prints, in order: the UTC start time, the host, the git HEAD and the worktree state,
the toolchain, a full capped build, `check_sorries.sh`, `check_axioms.sh`,
`check_root_imports.py`, `check_theorems_sigs.py`, and the UTC end time. Every dated or state dependent line of an
audit pack comes from this script, so no build result is hand written.

`AUDIT_PACK_JOBS` sets the core count of `lake-build-capped.sh` (default 12). Set it to the
thread budget of the session. `--no-build` skips the build; use it for a dry run only,
because the two gates then read stale oleans.

The script always exits 0. Read the recorded exit codes: a non-zero code means the pack must
not claim that the library builds.

## check_theorems_sigs.py - the third gate

```sh
python3 scripts/check_theorems_sigs.py docs/THEOREMS.md
```

It takes every Lean signature quoted in a document and looks for it in the tree, with
whitespace normalized and docstrings and comments removed. A signature that elides a proof
term as `…` or `_` matches anything at that position, so it reports `NO VERBATIM MATCH` only
when even the elided form is absent. It prints the count of signatures, the count matched
verbatim, and one line per mismatch. `NO VERBATIM MATCH` is expected for a signature the
document abbreviates on purpose and is a defect only when that abbreviation is undocumented;
`NAME NOT IN TREE` is always a defect. State on 2026-09-03: 92 signatures, 80 verbatim, 12
documented abbreviations, exit 0 (91, 79, 12 before the L5 block of that day; 71, 59, 12
before the 16 blocks added on 2026-09-02, 87, 75, 12 before the 4 blocks of the L2 and L3
additions). `docs/THEOREMS.md` states which of its signatures are abbreviated and why.

## check_paper_edits.py - the sixth gate, the paper edits as exact replacements

```sh
python3 scripts/check_paper_edits.py                 # the checks
python3 scripts/check_paper_edits.py --apply FILE    # also write the amended paper to FILE
```

In this public copy the paper snapshot and `notes/paper_edits.md` are absent; the script
prints `SKIP` and exits 0. The gate ran on the development tree (`docs/AUDIT_DOC.md`,
Appendix A, and the build of record in `docs/TECHNICAL.md`).

The paper snapshot `main_paper.tex` is never edited. `notes/paper_edits.md` records each
proposed edit as a subsection `### Exact replacement Ek.j: ...` with two fenced `latex`
blocks, the old text and the new text. `docs/THEOREMS.md` quotes each amended statement
under a lead-in `**Paper, as amended (Ek).**` followed by one fenced `latex` block. The
script reads no Lean. It checks that every old text occurs exactly once in the snapshot and
that no two overlap, that every new text has balanced braces and matched `\begin` and
`\end` pairs, and that every amended block of `THEOREMS.md` is verbatim in the amended paper
(the snapshot with every replacement applied, held in memory) and absent from the
unamended snapshot. It prints one line per hunk, in paper order, and one summary line; it
exits 1 on the first failure. `--apply` writes the amended paper only to a path outside the
repository or containing `scratch`; the amended paper is not a repository file, by the
user's decision of 2026-09-05. It does not compile the LaTeX (the class file and
`tdefs.tex` are not in the repository), so a macro that the new text uses must already
occur in the paper; the review of each hunk checks this by hand. State on 2026-09-05: 24
replacements, 6 amended blocks, exit 0, 0.1 s.

## check_core_imports.py - the eighth gate, the import boundary of the Gaussian RMT core

```sh
python3 scripts/check_core_imports.py              # gate, under 1 s
python3 scripts/check_core_imports.py --list       # every core module with its line count
python3 scripts/check_core_imports.py --boundary   # the paper-specific imports of the Gaussian
                                                   # discharge families outside the core
```

The rank-one Gaussian random matrix theory (`Defs`, `Spectral`, `RMT`, `RMT/*.lean` without
`RMT/Het/`, `Prob/`, six files of `LinAlg/`, `Vendor/COLT83/`; entry theorem
`singleTableLaw_of_gaussian` in `RMT/Full.lean`) is meant to be reusable on its own
(`README.md`, "The Gaussian RMT core"). The script fixes that module set in its source, reads
the `import StackedSVD...` lines of every core file, and fails when one of them imports a
module outside the set, or when the entry theorem is missing. Exit codes: 0 pass, 1 a
boundary violation or the missing entry theorem, 2 environment error. `--boundary` is
informational: it lists, for `RMT/Het/`, `RankR/RMT/`, `RankR/Het/` and
`RankR/SingleWeight/Het/`, the imports that are neither in the core nor in one of those four
families, which is the list the F28 refactor has to move.

## check_kernel.sh - the seventh gate, the kernel replay

```sh
scripts/check_kernel.sh                                    # every module of the package, about 10 min
scripts/check_kernel.sh StackedSVD.Defs StackedSVD.RankR   # only these module-name prefixes
KERNEL_CHECK_JOBS=2 scripts/check_kernel.sh                # a smaller worker pool
```

`lake build` sends each declaration through the kernel once, at elaboration time. The
toolchain's `leanchecker` (Lean v4.28 and later; the archived `lean4checker` before that)
reads the compiled `.olean` of a module, imports the oleans of its imports, and adds every
constant of the module to that environment through the kernel a second time
(`Lean.Environment.replay`), with no elaborator, no tactic and no `unsafe` code in the way.
It catches a declaration that reached the olean without a kernel check (environment hacking:
`addDeclCore (doCheck := false)`, `debug.skipKernelTC`, a metaprogram that edits the
environment) and a corrupt olean. It trusts the imports as loaded, so the gate replays every
module of ours and trusts Mathlib and StatsMLlib as compiled. It is not an external verifier:
the kernel that re-checks is the one of the same toolchain. Gate 2 does not cover this case:
`#print axioms` of a bogus theorem whose value is a constant prints no axiom at all.

Negative test (2026-09-07, scratch package outside the repository, same toolchain): a module
with the `AddFalse` example of the Lean test suite, `false : False` added with
`addDeclCore (doCheck := false)`, passes `lake build` (`Build completed successfully (5 jobs)`)
and fails the gate: `leanchecker found a problem in Kneg.AddFalse`, `(kernel) declaration
type mismatch, 'false' has type Prop but it is expected to have type False`, exit 1. A sound
module next to it replays with exit 0.

On the development server the gate runs leanchecker on the root-disk olean copy of
`lean-local.sh` (`LEAN_LOCAL_LIBS`, default `/tmp/$USER-lean-libs`; import 7 s per module
instead of 5 to 10 min over NFS) after it has compared every `.olean` of the
package in the copy with the build output byte for byte (`cmp`); a difference means the copy
is stale, exit 2, run `lean-local.sh sync`. Without the copy it uses the search path of `lake
env`. leanchecker starts one task per matched module at once and each task holds about 5 GB,
so on the development server the script caps the worker pool with the nprocs shim of
`lake-build-capped.sh` (`KERNEL_CHECK_JOBS`, default 3); without the cap the 182 tasks of the
package reached 231 GB on 2026-09-07. The shim exists only on that server. Elsewhere the
script runs `scripts/KernelReplay.lean` instead: a copy of the toolchain's `LeanChecker.lean`
(Apache 2.0, the original header kept) whose only change of substance is a worker pool
bounded by `-j`, run as `lean -j N --run scripts/KernelReplay.lean -v -j N <targets>` (the
first `-j` sizes the thread pool of the Lean runtime, which otherwise follows the CPU count
and reserves more address space than a 40 GB `ulimit -v` allows; the second sizes the pool
of replays). `KERNEL_CHECK_JOBS` sets both. Measured on 2026-09-10 in a fresh clone of the
public copy: 230 modules in 737 s with 3 workers, peak 20 GB. That count is 230, not 231,
because a fresh `lake build` writes no olean for `Vendor/COLT83/Axioms.lean`, the
`#print axioms` check file that nothing imports; the development tree has one from an
explicit build. The verdict line counts the `replaying` lines against the oleans of the
package: `check_kernel: OK - N modules replayed through the kernel, 0 problems; every one of
the N oleans of the package (S s, 3 workers)`.

## The numeric checks

These six scripts are not gates. Each one recomputes a closed form or a limit outside Lean
and asserts it, so a Lean statement that quotes a constant has an independent witness. They
need numpy. Each prints its seed. Each one
caps the BLAS threads at 4 (`OMP_NUM_THREADS` and its siblings, set before numpy loads,
since 2026-09-10 evening: an uncapped run of `check_nongaussian_forms.py` took 60 cores);
export those variables first to choose another count.

| Script | What it checks | Runtime |
|---|---|---|
| `check_edge_count.py` | task U7 of `notes/archive/rankr_plan_A.md`: the edge and the outlier count of the rank-`r` subspace law. Seed 2026090101. **Exits 1 today**: its C2 bar is an exact empirical rate of 1.000, which no finite `d` reaches near the edge. Read the cross-`d` trend block, and `notes/archive/audit_numeric_edge_count_2026-09-02.md` for the 32 rows and why none of them contradicts a claim. | 368 s |
| `check_port_stacksvd_subspace.py` | `prop_stacksvd_subspace_general` (`notes/archive/rankr_plan_B.md` section 6). Rows P1 to P5 are exact identities and must hit machine precision. | seconds |
| `check_rayleigh_bound.py` | the Rayleigh quotient bound uniform in the weights, rank 1 and rank `r` (paper findings E1 and E2). | seconds |
| `check_remark_examples.py` | the four examples of `remark:stack_outperform_svd` and `remark:svd_outperform_stack` (`main_paper.tex` 566 to 620). Deterministic, no seed. | seconds |
| `check_singleweight.py` | the single-weight stackSVD limit of Appendix D (`prop:gen_rank_stacksvd_singleweight`, paper finding E11): the closed form against Monte Carlo, and the value at a tied root. Seed 20260905. Prints its rows; no PASS bar. | about 2 minutes on 2 cores |
| `check_stacksvd_orthogonality.py` | the asymptotic orthogonality of the rank-`r` svdstack columns (`thm_rank_r_svdstack_offdiag_of_sep`, follow-up F2). Seed 20260905. No PASS bar. | seconds |
| `weighted_edge_levels.py` | the edge levels of the weighted noise block (`bHet`, `bSF`, the Weyl and operator-norm bounds of `notes/WEIGHTED_GENERAL_SCOPE.md` section 4): the share of supercritical draws whose outlier lies above each level. Seed 20260910. No PASS bar. | about 20 s |

The scripts `check_singleweight_regimes.py` and `check_nongaussian_forms.py` (the three weight
regimes of F18b; the non-Gaussian resolvent forms, which reads the R and Python reference code
of the paper from a checkout of the `stackedSVD` repository, `STACKEDSVD_REPO` or the parent
of this tree) are the same kind of witness. `scripts/numeric/` (the top-level
`scratchpad/` until 2026-09-09) holds the exploration scripts behind two design decisions,
with its own README; they are evidence, not checks.

## scripts/stage3_count/ - numeric evidence for the Furedi-Komlos count

Six scripts, kept as evidence for `cellCard_mul_le_pos`, the Furedi-Komlos count of Stage 3
(proved 2026-09-10 evening, `lean/StackedSVD/RMT/General/Edge/CountBound.lean`, moved from
`Count.lean`). Exact integer enumeration, no randomness, runtimes under 1 s to 17 s. Not a
gate; see `scripts/stage3_count/README.md` for what each script checks and the result.
Three more scripts (`x8_enum.py`, `x8_charge.py`, `x8_cells8.py`) back the X8 plan that
proved the same lemma, `notes/x8_plan.md` (a six-part code of a walk and a Narayana lower
bound on the tree cell). `scripts/stage3_count/x8_scratch/` (development repository only, not
part of this copy) holds the scratch Lean files and reports of that plan's units, produced
2026-09-10, not part of the Lake build.

## check_paths.py - every path the documents name exists

Not a gate of the build of record; a check on the documents (added 2026-09-10 with the
repository cleanup, `notes/REPO_CLEANUP_PLAN.md`). It reads `README.md`, `CLAUDE.md`,
`docs/*.md`, `scripts/README.md`, `scripts/numeric/README.md` and `lean/README.md`,
takes every backtick-quoted token that looks like a path (a slash, or a file extension), and
resolves it against the repository root, the Lean package, the Lean sources, `docs/`,
`notes/`, `notes/archive/`, `scripts/`, and finally as a path suffix or a bare file name
anywhere in the tracked tree. A mention of a deleted file survives when the same line says
so (`deleted`, `removed`, `retired`, `not kept`, `until`, `was`, `were`, a date). Exit 0
when every path resolves (the count of record is in `docs/TECHNICAL.md`), 1 with the list of
the missing ones.

```sh
python3 scripts/check_paths.py            # gate form
python3 scripts/check_paths.py --list     # every token and where it resolved
```

Release mode (2026-09-10): the public copy of the tree ships without `notes/`, the paper
snapshot and a few private documents. When `scripts/release_withheld.txt` exists (the release
script writes it; the development tree has none, and the file is not tracked), a token that does not resolve in the tree
but resolves against that list counts as withheld, not as missing, and the summary line
reports the count. Outside a git checkout the script walks the tree instead of `git ls-files`.

## release/make_release.py - the public copy

```sh
python3 scripts/release/make_release.py OUTDIR [--tar]   # writes OUTDIR/lean_formalization/
```

Exports the public copy (the folder `lean_formalization` of the paper's code repository,
`https://github.com/phillipnicol/stackedSVD`; the user's decision of 2026-09-10, the plan is
`notes/RELEASE_PLAN.md`). It copies the tracked files that ship (the Lean project, the reader
documents, the gate and numeric scripts, `README.md`, `GETTING_STARTED.md`, `LICENSE`,
`CITATION.cff`), applies a fixed list of exact text replacements (each must match once: the
clone lines, the server paths, the sentences that name a withheld file), adds an "About this
copy" section to `README.md`, writes the list of the withheld paths (`release_withheld.txt`
under `scripts/`, generated, not tracked), and then checks the
copy: every unmodified file byte-identical, no server path or session link in the reader
documents, `check_paths.py` inside the copy with 0 missing. Not shipped: `notes/` (which since
2026-09-11 holds the paper snapshot and the two paper-facing documents too), `CLAUDE.md`,
`scripts/stage3_count/x8_scratch/`, `scripts/make_audit_packet_rank1.py` and the folder
`scripts/release/` itself. The docstring lists the rules. `release/make_branch.sh` puts the
exported folder on a branch of a clone of the paper's code repository (no push), and
`release/README.md` gives the commands from a clone.

## AuditSignatures.lean

Not part of the build graph. It prints the elaborated signature of every headline theorem,
with the ambient variables and instances that a `section` hides in the source, and
`#print axioms` for each. It answers action items 14 and 17 of the external audit of
2026-08-30. Run it with `lake env lean scripts/AuditSignatures.lean` from `lean/`.
