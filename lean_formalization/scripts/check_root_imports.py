#!/usr/bin/env python3
"""Gate: every module of the package is reachable from the root import file.

`scripts/AxiomAudit.lean` audits what `import StackedSVD` pulls in. A module that
`lean/StackedSVD.lean` does not reach, directly or transitively, is compiled by
`lake build` and then never audited, so a `sorry` or an `axiom` in it passes both gates
unseen. This script closes that hole.

It reads the import lines of every `lean/StackedSVD/**/*.lean` file, including the
`public import` and `private import` forms of the Lean 4.33 module system, computes the
transitive closure from the root file, and fails when a module of ours is outside it.

Vendored code (`StackedSVD/Vendor/`) is reported but not gated: `AxiomAudit.lean` skips it
unless `AXIOM_AUDIT_VENDOR=1`, and `Vendor/COLT83/Axioms.lean` is the vendor's own
`#print axioms` driver, which nothing imports on purpose.

Usage (from anywhere):
  python3 scripts/check_root_imports.py            gate
  python3 scripts/check_root_imports.py --list     print the closure size and every module

Exit codes: 0 pass, 1 a module of ours is unreachable, 2 environment error.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PROJ = os.path.join(ROOT, "lean")
PKG = os.path.join(PROJ, "StackedSVD")
ROOT_FILE = os.path.join(PROJ, "StackedSVD.lean")

IMPORT = re.compile(r"^(?:public\s+|private\s+|meta\s+)*import\s+(StackedSVD\S*)")


def imports_of(path):
    out = set()
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            m = IMPORT.match(line.strip())
            if m:
                out.add(m.group(1))
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "gate"
    if mode not in ("gate", "--list"):
        print("check_root_imports.py: unknown option %r" % mode, file=sys.stderr)
        return 2
    if not os.path.isdir(PKG) or not os.path.isfile(ROOT_FILE):
        print("check_root_imports.py: no package at %s" % PKG, file=sys.stderr)
        return 2

    mods = {}
    for dirpath, _dirnames, filenames in os.walk(PKG):
        for f in sorted(filenames):
            if f.endswith(".lean"):
                path = os.path.join(dirpath, f)
                name = os.path.relpath(path, PROJ)[:-5].replace(os.sep, ".")
                mods[name] = imports_of(path)

    seen, stack = set(), list(imports_of(ROOT_FILE))
    while stack:
        cur = stack.pop()
        if cur in seen:
            continue
        seen.add(cur)
        stack.extend(mods.get(cur, set()) - seen)

    ours = {m for m in mods if not m.startswith("StackedSVD.Vendor.")}
    vendor = set(mods) - ours
    miss_ours = sorted(ours - seen)
    miss_vendor = sorted(vendor - seen)

    if mode == "--list":
        for m in sorted(mods):
            print("%-4s %s" % ("ok" if m in seen else "MISS", m))

    print("check_root_imports: %d modules (%d ours, %d vendored); %d reachable from "
          "lean/StackedSVD.lean." % (len(mods), len(ours), len(vendor),
                                           len(seen & set(mods))))
    if miss_vendor:
        print("  vendored and unreachable (not gated, see the module docstring):")
        for m in miss_vendor:
            print("    %s" % m)
    if not miss_ours:
        print("check_root_imports: OK - every module of the package is audited.")
        return 0
    print("check_root_imports: FAIL - %d module(s) of ours are outside the closure, so "
          "AxiomAudit.lean never sees them:" % len(miss_ours))
    for m in miss_ours:
        print("    %s" % m)
    print("  Add each to lean/StackedSVD.lean, or import it from a module that is "
          "already in the closure.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
