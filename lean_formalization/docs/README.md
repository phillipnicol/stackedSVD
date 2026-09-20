# docs/: the detailed documents

Four documents. `README.md` at the repository root is
the entry point for a reader; this folder holds the detail behind it. Every document here is
about the formalization; none of them changes the paper.

| File | What it is |
|---|---|
| `TECHNICAL.md` | the detailed map: the notation table (paper symbol to Lean name), the main-results table, how the random matrix theory input enters, verification, the layout of the tree file by file, and the build of record with its timings. It was the top-level `README.md` until 2026-09-03. |
| `THEOREMS.md` | the self-contained auditor annex. For every result: the paper statement, the Lean statement, the hypotheses in words, the modeling choices, the scope, the audit trail. Gate 4 (`scripts/check_theorems_sigs.py`) checks that the quoted signatures match the source. |
| `AUDIT_DOC.md` | the self-contained audit entry point. Read it without the repository open and you learn what is claimed, what is assumed, what is not proved, and how to check every claim. Appendix A holds the transcript of the build and the eight gates. |
| `SORRIES.md` | the ledger of open proofs. Every `sorry` under `lean/StackedSVD/` needs a row here. The table is empty, and gate 1 (`scripts/check_sorries.sh`) keeps it so. |

`GETTING_STARTED.md` (install Lean, build the project, check one theorem; by Phillip Nicol)
lived in this folder until 2026-09-10. It is now at the repository root, where a first-time
reader finds it next to `README.md`.

`paper/main_paper.tex` (the read-only paper snapshot), `paper_edits.md` (12 findings about
the paper, proposals only) and `NOTE_FOR_COAUTHORS.md` (one page for a coauthor) lived in this
folder until 2026-09-11. They are now in `notes/`, since they are private and do not ship in
the public release.
