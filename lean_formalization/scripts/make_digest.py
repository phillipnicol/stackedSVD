#!/usr/bin/env python3
"""Signature digest of Lean files: module docs, docstrings, and declaration headers, proofs cut.

Usage: python3 scripts/make_digest.py [--exported [--tree DIR]] <file.lean> [...] > digest.md
Each declaration header is printed from its keyword line to the first `:=` (the `:=` and the
proof are dropped); a `structure` or `class` is printed to the next blank line, so its fields
show. Lines of a `section`, `namespace`, `variable`, `open`, `omit` and `set_option` are kept,
since they change the meaning of the headers below them. Nothing is elaborated: this is text.

`--exported`: a declaration whose short name occurs in another file of the tree (default
`lean/StackedSVD`, or `--tree DIR`) keeps its docstring and header and gets a
`-- used by: k files (...)` tag; every other declaration is listed by name only, with the
number of its uses inside its own file. Module docs (`/-! ... -/`) are always kept.
"""
import os
import re
import sys
from collections import Counter

KEYWORDS = re.compile(
    r"^(?:@\[[^\]]*\]\s*)?(?:private |protected |noncomputable |nonrec )*"
    r"(theorem|lemma|def|abbrev|instance|structure|class|inductive|opaque|axiom)\s+"
    r"([A-Za-z_][\w.']*)?")
KEEP = re.compile(r"^(section|end|namespace|variable|open|omit|set_option|universe|import)\b")
IDENT = re.compile(r"[A-Za-z_][\w']*")


def token_sets(tree):
    toks = {}
    for r, _, fs in os.walk(tree):
        for f in fs:
            if f.endswith(".lean"):
                p = os.path.join(r, f)
                toks[p] = Counter(IDENT.findall(open(p).read()))
    return toks


def users(name, path, toks):
    short = name.split(".")[-1]
    me = os.path.abspath(path)
    return sorted(os.path.relpath(p, "lean/StackedSVD").removesuffix(".lean")
                  for p in toks if os.path.abspath(p) != me and short in toks[p])


def digest(path, toks=None):
    out = [f"\n## {path} ({sum(1 for _ in open(path))} lines)\n"]
    lines = open(path).read().split("\n")
    i, n = 0, len(lines)
    pending = []  # docstring lines that belong to the next declaration
    while i < n:
        line = lines[i]
        s = line.strip()
        if s.startswith("/-!"):
            while i < n:
                out.append(lines[i])
                if "-/" in lines[i]:
                    i += 1
                    break
                i += 1
            continue
        if s.startswith("/--"):
            pending = []
            while i < n:
                pending.append(lines[i])
                if "-/" in lines[i]:
                    i += 1
                    break
                i += 1
            continue
        if KEEP.match(s):
            out.append(line)
            i += 1
            continue
        m = KEYWORDS.match(s)
        if m:
            kind, name = m.group(1), m.group(2) or "?"
            head = []
            if kind in ("structure", "class", "inductive"):
                while i < n and lines[i].strip() != "":
                    head.append(lines[i])
                    i += 1
            else:
                while i < n:
                    cur = lines[i]
                    if ":=" in cur:
                        h = cur[: cur.index(":=")].rstrip()
                        if h.strip():
                            head.append(h)
                        i += 1
                        break
                    head.append(cur)
                    i += 1
            if toks is None:
                out.extend(pending)
                out.extend(head)
                out.append("")
            else:
                us = users(name, path, toks)
                if us or kind in ("structure", "class", "inductive"):
                    out.extend(pending)
                    out.extend(head)
                    shown = ", ".join(us[:6]) + (", ..." if len(us) > 6 else "")
                    out.append(f"-- used by: {len(us)} files" + (f" ({shown})" if us else ""))
                    out.append("")
                else:
                    own = toks[os.path.abspath(path)][name.split(".")[-1]] - 1 if os.path.abspath(path) in toks else -1
                    out.append(f"-- helper {kind} {name}: {own} uses in this file only")
            pending = []
            continue
        if s == "" :
            pending = []  # a blank line detaches a docstring from what follows
        i += 1
    return "\n".join(out)


if __name__ == "__main__":
    args = sys.argv[1:]
    toks = None
    tree = "lean/StackedSVD"
    if "--tree" in args:
        k = args.index("--tree")
        tree = args[k + 1]
        del args[k:k + 2]
    if "--exported" in args:
        args.remove("--exported")
        toks = {os.path.abspath(p): c for p, c in token_sets(tree).items()}
    for p in args:
        sys.stdout.write(digest(p, toks))
