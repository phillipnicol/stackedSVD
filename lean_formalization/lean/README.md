# StackedSVD (the Lake package)

This directory is the Lean 4 package: `lakefile.toml` (Mathlib v4.33.0 and StatsMLlib,
both pinned), `lean-toolchain`, the root module `StackedSVD.lean`, and the source tree
`StackedSVD/`. Every `lake` command runs from here.

The documentation lives one level up: `../README.md` is the entry point,
`../GETTING_STARTED.md` the install-and-check walkthrough, `../docs/TECHNICAL.md` the
map of the source tree, and `../docs/THEOREMS.md` the statement of every theorem next to the
paper's. `check_weighted.lean` is the four-line example of the walkthrough (`lake env lean
check_weighted.lean` prints the type of one theorem).
