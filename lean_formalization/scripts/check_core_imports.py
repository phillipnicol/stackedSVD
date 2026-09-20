#!/usr/bin/env python3
"""Gate: the Gaussian RMT core, and the wider set EXT, import only modules of their own kind.

The rank-one Gaussian random matrix theory of this project (the Marchenko-Pastur resolvent
estimates, the Gaussian concentration and symmetry arguments, and their conclusion
`SpikedModel.singleTableLaw_of_gaussian`, the BBP-type limits of the Gaussian rank-one
spiked model) is meant to be reusable by a reader who does not care about stacked SVD. This
script fixes two sets of modules, the core and the wider set EXT, and fails when a module
imports outside its own set, when a listed EXT module is missing, or when an entry theorem
is gone.

The core (2026-09-08): `Defs`, `Spectral`, `RMT` (the hypothesis structure `SingleTableLaw`),
every `RMT/*.lean` except `RMT/Het/`, every `Prob/*.lean`, six files of `LinAlg/` (Weyl and
Davis-Kahan for the top projector and the top-r spectral projector, Ky Fan, padded spectra,
top-r frames, the reindexing of the spectral projector) and the vendored `Vendor/COLT83/`
Gaussian files. `LinAlg/SpecIdx*.lean` and `LinAlg/SpecWindow.lean` are outside the core,
since they import `SVDStack/Defs.lean` for `vMax`; the heteroscedastic (`RMT/Het/`) and
rank-r (`RankR/RMT/`, `RankR/Het/`, `RankR/SingleWeight/Het/`) Gaussian discharges are
outside the core, since they import the paper's definition files (`StackSVDWeighted.lean`,
`RankR/General.lean`) for the structures they discharge. All of these join the wider set
EXT, described next.

EXT, the extended set (2026-09-08 widening, item F28 of `notes/FOLLOWUP_LIST.md`): the
core, plus every module of the four families above except their three `Sup` files, plus 21
named modules that hold the paper's model definitions and the Layer 0 stacking lemmas
(`SVDStack/Defs`, `Scalars`, `StackSVD`, `RankR/Defs`, and so on; `EXT_EXPLICIT` below lists
all 21). EXT is import-closed: every import line of an EXT module names another EXT module.
No declaration inside EXT takes one of the 8 law structures (`SingleTableLaw`, `HeteroLaw`,
`ThetaEstLaw`, `TableLawR`, `SubspaceLaw`, `SubspaceLawG`, `HeteroLawR`, `SingleWeightLaw`)
as a hypothesis, so EXT stays free of Layer 1 consumers. This script checks only the import
closure. Gate 5 (`scripts/check_layering.sh`) checks the consumer-free property; this
script does not.

The three `Sup` files (`RMT/Het/Sup.lean`, `RankR/Het/Sup.lean`,
`RankR/SingleWeight/Het/Sup.lean`) stay outside EXT. Each holds Layer 1 corollaries (for
example `thm_stacksvd_weighted_gaussian`, `thm_rank_r_stacksvd_gaussian`) that take a law
structure as a hypothesis by design, so they cannot join a consumer-free set. `--boundary`
prints their imports, and the imports of every other non-EXT module that imports an EXT
module.

Usage (from anywhere):
  python3 scripts/check_core_imports.py                gate on the real package
  python3 scripts/check_core_imports.py --pkg DIR       gate on a package copy at DIR
  python3 scripts/check_core_imports.py --list          print the core, then the EXT
                                                        modules that are not core, each
                                                        with its line count
  python3 scripts/check_core_imports.py --boundary      print the imports of the 3 `Sup`
                                                        files, then of every other non-EXT
                                                        module that imports an EXT module

Exit codes: 0 pass, 1 a core or EXT module imports outside its own set, a listed EXT module
is missing, or an entry theorem is missing, 2 environment error.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DEFAULT_PKG = os.path.join(ROOT, "lean", "StackedSVD")

IMPORT = re.compile(r"^(?:public\s+|private\s+|meta\s+)*import\s+(StackedSVD\S*)")
ENTRY = re.compile(r"^theorem singleTableLaw_of_gaussian\b")

CORE_SINGLE = {"StackedSVD.Defs", "StackedSVD.Spectral", "StackedSVD.RMT"}
CORE_PREFIX = ("StackedSVD.RMT.", "StackedSVD.Prob.", "StackedSVD.Vendor.")
CORE_EXCLUDE_PREFIX = ("StackedSVD.RMT.Het",)
CORE_LINALG = {
    "StackedSVD.LinAlg." + x
    for x in ("Eigen", "Frame", "KyFan", "SpecInvReindex", "SpecProjPerturb", "TopProjPerturb")
}

# The four Gaussian discharge families that join EXT, each minus its own Sup corollary file.
EXT_FAMILY_PREFIX = (
    "StackedSVD.RMT.Het.",
    "StackedSVD.RankR.RMT.",
    "StackedSVD.RankR.Het.",
    "StackedSVD.RankR.SingleWeight.Het.",
)
SUP_MODULES = (
    "StackedSVD.RMT.Het.Sup",
    "StackedSVD.RankR.Het.Sup",
    "StackedSVD.RankR.SingleWeight.Het.Sup",
)

# The 21 named modules of notes/F28_PLAN.md section 2: the paper's model definitions and
# the Layer 0 stacking lemmas that EXT needs, that no family prefix rule already reaches.
EXT_EXPLICIT = (
    "StackedSVD.LinAlg.SpecIdx",
    "StackedSVD.LinAlg.SpecIdxMeas",
    "StackedSVD.LinAlg.SpecIdxPerturb",
    "StackedSVD.LinAlg.SpecWindow",
    "StackedSVD.SVDStack.Defs",
    "StackedSVD.SVDStack.DelocDir",
    "StackedSVD.Scalars",
    "StackedSVD.StackSVD",
    "StackedSVD.StackSVDWeighted",
    "StackedSVD.RankR.Defs",
    "StackedSVD.RankR.General",
    "StackedSVD.RankR.Flatten",
    "StackedSVD.RankR.Weighted",
    "StackedSVD.RankR.Subspace",
    "StackedSVD.RankR.SubspaceG",
    "StackedSVD.RankR.SubspaceGStack",
    "StackedSVD.RankR.Aligned",
    "StackedSVD.RankR.GramR",
    "StackedSVD.RankR.StackGamma",
    "StackedSVD.RankR.SingleWeight.Scalars",
    "StackedSVD.RankR.SingleWeight.Defs",
)
assert len(EXT_EXPLICIT) == 21, len(EXT_EXPLICIT)

# The four entry theorems of the Gaussian discharges (grep of 2026-09-08). Three of their
# files sit outside EXT (they are Sup corollary files); the check here is only that the
# named theorem exists in the named file, the same as the core's singleTableLaw_of_gaussian.
ENTRY_THEOREMS = (
    ("heteroLaw_of_gaussian", ("RMT", "Het", "Sup.lean")),
    ("tableLawR_of_gaussian_rk", ("RankR", "RMT", "TableLawGaussian.lean")),
    ("heteroLawR_of_gaussian", ("RankR", "Het", "Sup.lean")),
    ("singleWeightLaw_of_gaussian", ("RankR", "SingleWeight", "Het", "Sup.lean")),
)


def is_core(mod):
    if mod in CORE_SINGLE or mod in CORE_LINALG:
        return True
    if mod.startswith(CORE_EXCLUDE_PREFIX):
        return False
    return mod.startswith(CORE_PREFIX)


def is_ext(mod):
    if is_core(mod):
        return True
    if mod in EXT_EXPLICIT:
        return True
    return mod.startswith(EXT_FAMILY_PREFIX) and mod not in SUP_MODULES


def modules(pkg):
    """Every module of the package at pkg: name -> (path, imports, line count)."""
    out = {}
    for dp, _, fs in os.walk(pkg):
        for f in sorted(fs):
            if not f.endswith(".lean"):
                continue
            path = os.path.join(dp, f)
            rel = os.path.relpath(path, pkg)[:-5]
            name = "StackedSVD." + rel.replace(os.sep, ".")
            imps = []
            n = 0
            with open(path, encoding="utf-8") as fh:
                for line in fh:
                    n += 1
                    m = IMPORT.match(line.strip())
                    if m:
                        imps.append(m.group(1))
            out[name] = (path, imps, n)
    return out


def counts(mods, names):
    """(modules present, their lines, ours among them, ours' lines) for a set of names."""
    present = [m for m in names if m in mods]
    ours = [m for m in present if not m.startswith("StackedSVD.Vendor.")]
    lines = sum(mods[m][2] for m in present)
    lines_ours = sum(mods[m][2] for m in ours)
    return len(present), lines, len(ours), lines_ours


def main(argv):
    pkg = DEFAULT_PKG
    if "--pkg" in argv:
        i = argv.index("--pkg")
        if i + 1 >= len(argv):
            print("check_core_imports: --pkg needs a directory argument", file=sys.stderr)
            return 2
        pkg = argv[i + 1]
    if not os.path.isdir(pkg):
        print(f"check_core_imports: package directory not found: {pkg}", file=sys.stderr)
        return 2

    mods = modules(pkg)
    core = sorted(m for m in mods if is_core(m))
    ext = sorted(m for m in mods if is_ext(m))
    if not core:
        print("check_core_imports: no core module found", file=sys.stderr)
        return 2

    n_core, lines_core, ours_core, lines_core_ours = counts(mods, core)
    n_ext, lines_ext, ours_ext, lines_ext_ours = counts(mods, ext)

    if "--list" in argv:
        for m in core:
            print(f"{mods[m][2]:6d}  {m}")
        extra = [m for m in ext if m not in core]
        print(f"-- EXT extra, not core ({len(extra)} modules) --")
        for m in extra:
            print(f"{mods[m][2]:6d}  {m}")

    if "--boundary" in argv:
        boundary = set(SUP_MODULES)
        for m in mods:
            if not is_ext(m) and any(is_ext(i) for i in mods[m][1]):
                boundary.add(m)
        print(f"boundary ({len(boundary)} modules): the 3 Sup files, plus every other")
        print("non-EXT module that imports an EXT module")
        for m in sorted(boundary):
            tag = " (Sup)" if m in SUP_MODULES else ""
            print(f"  {m}{tag} imports:")
            for i in sorted(mods.get(m, ("", [], 0))[1]):
                print(f"    {i}")

    ok = True

    bad_core = [(m, i) for m in core for i in mods[m][1] if not is_core(i)]
    bad_ext = [(m, i) for m in ext for i in mods[m][1] if not is_ext(i)]
    missing_import = sorted({i for m in ext for i in mods[m][1] if i not in mods})
    missing_listed = [m for m in EXT_EXPLICIT if m not in mods]

    entry_ok = False
    try:
        with open(os.path.join(pkg, "RMT", "Full.lean"), encoding="utf-8") as fh:
            entry_ok = any(ENTRY.match(line) for line in fh)
    except OSError:
        pass

    print(f"core: {n_core} modules, {lines_core} lines ({ours_core} ours, "
          f"{lines_core_ours} lines; {n_core - ours_core} vendored)")
    print(f"EXT: {n_ext} modules, {lines_ext} lines ({ours_ext} ours, "
          f"{lines_ext_ours} lines; {n_ext - ours_ext} vendored)")

    for m, i in bad_core:
        print(f"FAIL: core module {m} imports non-core module {i}")
        ok = False
    for m, i in bad_ext:
        print(f"FAIL: EXT module {m} imports non-EXT module {i}")
        ok = False
    for i in missing_import:
        print(f"FAIL: import of unknown module {i}")
        ok = False
    for m in missing_listed:
        print(f"FAIL: listed EXT module {m} not found under {pkg}")
        ok = False

    if entry_ok:
        print("entry theorem singleTableLaw_of_gaussian: present in RMT/Full.lean")
    else:
        print("FAIL: theorem singleTableLaw_of_gaussian not found in RMT/Full.lean")
        ok = False

    for name, parts in ENTRY_THEOREMS:
        path = os.path.join(pkg, *parts)
        pat = re.compile(r"^theorem " + re.escape(name) + r"\b")
        found = False
        try:
            with open(path, encoding="utf-8") as fh:
                found = any(pat.match(line) for line in fh)
        except OSError:
            found = False
        rel = "/".join(parts)
        if found:
            print(f"entry theorem {name}: present in {rel}")
        else:
            print(f"FAIL: theorem {name} not found in {rel}")
            ok = False

    print("check_core_imports: PASS" if ok else "check_core_imports: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
