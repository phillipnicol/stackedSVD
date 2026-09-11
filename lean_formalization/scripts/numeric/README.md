# scripts/numeric/ - numeric exploration, not part of any gate

These scripts are the working numerics behind two design decisions. They are kept as evidence
of what was checked, they are not gates, and nothing in the Lean tree depends on them. The
gates and the maintained numeric checks are in `scripts/` (see `scripts/README.md`). Until
2026-09-09 this directory was the top-level `scratchpad/`; the notes written before that
date name the old path.

Until 2026-09-10, 15 docstrings in 11 Lean files cited numeric checks as
`scratchpad/<unit>/<script>.py`; those were session scratch scripts that were not kept, and
the docstrings now say so (the one survivor, `check_trackG.py`, is tracked at
`notes/archive/trackG_specs/check_trackG.py` and cited by that path). The printed numbers a
docstring reports are recorded in the review note of that unit under `notes/archive/`. The 14
files here are the ones that were kept.

| Path | What it explored |
|---|---|
| `audit_p1_rayleigh.py` | the Rayleigh quotient bound uniform in the weights, at `M = 4`, `d = 60`, seed 20260901. The maintained version is `scripts/numeric/check_rayleigh_bound.py`. |
| `heteroedge_elementary/` | the heteroscedastic bulk edge: the scalars `zfun`, `phi`, `sLo`, `sStar`, `bHet`, `gammaTop` transcribed from `RMT/Het/MPhet.lean` and `StackSVDWeighted.lean`, and Monte Carlo sweeps against them (`mc.py`, `sweep.py`, `bigsweep.py`, `rev.py`, `split.py`, `verify_stair.py`, and the `report_*.py` readers). |

They need numpy. Each prints its seed. Written 2026-08-30 and 2026-09-01.
