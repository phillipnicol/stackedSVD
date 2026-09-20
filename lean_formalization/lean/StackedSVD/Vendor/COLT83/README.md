# Vendored COLT83 Gaussian files

Source: <https://github.com/RemyDegenne/colt-2026-83>, commit
`2d846025a86ee6760a805901c04237039f987457`, Apache-2.0 (see `LICENSE`).
Written for Mathlib v4.34.0-rc2; backported to the pinned Mathlib v4.33.0 on 2026-08-29.

19 modules, the transitive in-repo closure of `SteinReal`, `SteinIdentity`,
`GaussianInterpolation`, `SudakovFernique` and `BorellTIS`. The path under this directory
matches the upstream path under `COLT83/`.

Three edits per file, and nothing else:

1. One comment line under the license header that records the source commit.
2. `import COLT83.…` becomes `import StackedSVD.Vendor.COLT83.…`.
3. `set_option autoImplicit false` after the imports. Upstream sets this in its lakefile;
   this project's lakefile does not, so the option moves into each file.

No statement or proof changed. No Mathlib name needed a rename between v4.34.0-rc2 and
v4.33.0. `Axioms.lean` is new; it prints the axioms of the two results the project uses.
