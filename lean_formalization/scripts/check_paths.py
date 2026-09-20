#!/usr/bin/env python3
"""check_paths.py - every path named in the reader-facing documents exists.

Usage (from anywhere):
  python3 scripts/check_paths.py            gate: exit 0 if every path resolves, 1 if not
  python3 scripts/check_paths.py --list     print every path it checked and where it resolved

Scope: README.md, CLAUDE.md, GETTING_STARTED.md, docs/*.md, scripts/README.md,
scripts/numeric/README.md, lean/README.md, notes/README.md, notes/discoveries/README.md,
notes/audit_packets/*/README.md (the reader-facing documents; the other notes are history
and are not checked).

A path is a backtick-quoted token that contains a slash or ends in .md, .lean, .py, .sh,
.yml, .cff, .tex, .toml, .json or .txt, without a placeholder (`<...>`, `*`, `...`) and without
a space. Each is resolved in order against: the repository root; `lean/` (the lake project);
`lean/StackedSVD/` (the Lean sources, which the documents name by their path inside the
package); `docs/`; `notes/`; `notes/paper/`; `notes/archive/`; `scripts/`; `scripts/numeric/`;
the directory of the document. A
token that names a deleted file with a date or the word "deleted", "removed", "until", or
"was" within 120 characters before it (the end of the previous line included) or 60 after it is skipped, so history stays
mentionable. Absolute paths (`/tmp/...` and the like) and URLs are skipped.

Release mode: a public copy of this tree ships without `notes/`, the paper snapshot and a
few private documents (see `scripts/release_withheld.txt`, written by the release script and
absent in the development tree). When that list exists, a token that does not resolve in the
tree but resolves against the list counts as "withheld", not as missing: the citation is
correct, the file is not part of the release.

Exit codes: 0 every path resolves, 1 at least one does not, 2 usage error.
"""
import os, re, sys

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
docs = ["README.md", "CLAUDE.md", "GETTING_STARTED.md", "scripts/README.md",
        "scripts/numeric/README.md", "lean/README.md", "notes/README.md",
        "notes/SERVER_SETUP.md", "notes/discoveries/README.md", "notes/audit_packets/README.md",
        "notes/paper_edits.md", "notes/NOTE_FOR_COAUTHORS.md"]
docs = [d for d in docs if os.path.isfile(os.path.join(root, d))]
packets = os.path.join(root, "notes", "audit_packets")
docs += sorted("notes/audit_packets/" + d + "/README.md" for d in os.listdir(packets)
               if os.path.isfile(os.path.join(packets, d, "README.md"))) if os.path.isdir(packets) else []
external = {"env.sh", "lean-local.sh", "lake-build-capped.sh", "simulations.py", "theory_pred.R",
            "CONTEXT.md", "tdefs.tex", "lean4checker", "nprocs_shim.so",
            ".lake",  # `.lake/` is the build output, present only after a build
            "release_withheld.txt"}  # generated into the public copy by scripts/release/make_release.py
import subprocess
# every tracked or unignored file of the repository (`.lake/` is ignored, so Mathlib stays out)
all_files = subprocess.run(["git", "-C", root, "ls-files", "--cached", "--others", "--exclude-standard"],
                           capture_output=True, text=True).stdout.split("\n")
all_files = [f for f in all_files if f]
if not all_files:  # not a git checkout (a copied release folder): walk the tree, `.lake/` and `.git/` excluded
    for d, subdirs, files in os.walk(root):
        subdirs[:] = [x for x in subdirs if x not in (".lake", ".git", "__pycache__")]
        all_files += [os.path.relpath(os.path.join(d, f), root) for f in files]
withheld_list = os.path.join(here, "release_withheld.txt")
withheld = [l.strip() for l in open(withheld_list)] if os.path.isfile(withheld_list) else []
withheld = [w for w in withheld if w and not w.startswith("#")]

def path_set(files):
    """The files and every directory above them, for an existence test by name."""
    out = set()
    for f in files:
        f = os.path.normpath(f)
        out.add(f)
        while "/" in f:
            f = os.path.dirname(f); out.add(f)
    return out

tree_set = path_set(all_files)
withheld_set = path_set(withheld)
docs += sorted("docs/" + f for f in os.listdir(os.path.join(root, "docs")) if f.endswith(".md"))
bases = ["", "lean", "lean/StackedSVD", "docs", "notes", "notes/paper", "notes/archive", "scripts",
         "scripts/numeric"]
exts = (".md", ".lean", ".py", ".sh", ".yml", ".cff", ".tex", ".toml", ".json", ".txt")
skip_words = re.compile(r"(deleted|removed|retired|not kept|until|was |were |renamed|formerly|no longer|moved (out|from)|20\d\d-\d\d-\d\d)", re.I)
tok_re = re.compile(r"`([^`\n]+)`")

mode = sys.argv[1] if len(sys.argv) > 1 else "gate"
if mode not in ("gate", "--list"):
    print(__doc__); sys.exit(2)

def resolve(tok, docdir, files, fileset, on_disk=False):
    """Where `tok` resolves among `files` (a list of paths; `fileset` is their path_set), or None.
    `on_disk`: a path that exists in the tree but is ignored by git (`lean/.lake/`) also counts."""
    t = tok.rstrip("/").split(":")[0]  # `file.lean:123` names a line
    cands = [t] if "." in os.path.basename(t) else [t, t + ".lean"]
    for c in cands:
        for b in bases + [docdir]:
            p = os.path.normpath(os.path.join(b, c))
            if p in fileset or (on_disk and os.path.exists(os.path.join(root, p))):
                return p
    # a path suffix inside the Lean sources (`Het/Sup.lean` names one of three files; any hit counts)
    for c in cands:
        if c.endswith(".lean"):
            hits = [f for f in files if f.startswith("lean/StackedSVD/") and f.endswith("/" + c)]
            if hits:
                return hits[0] if len(hits) == 1 else hits[0] + f" (+{len(hits) - 1} more)"
    # a directory suffix (`Het/`) or a bare file name anywhere in the repository
    if tok.endswith("/"):
        hits = [f for f in files if ("/" + t + "/") in ("/" + f)]
        if hits:
            return os.path.dirname(hits[0]) + "/"
    if "/" not in t:
        hits = [f for f in files if os.path.basename(f) == t]
        if hits:
            return hits[0] if len(hits) == 1 else hits[0] + f" (+{len(hits) - 1} more)"
    return None

n_ok = n_bad = n_withheld = 0
bad = []
for doc in docs:
    path = os.path.join(root, doc)
    if not os.path.exists(path):
        continue
    docdir = os.path.dirname(doc)
    prev = ""
    for ln, line in enumerate(open(path, encoding="utf-8"), 1):
        for m in tok_re.finditer(line):
            tok = m.group(1).strip()
            if " " in tok or "<" in tok or "*" in tok or "..." in tok or "…" in tok:
                continue
            if tok.startswith(("/", "http", "~", "$", "-", "#", "@")) or tok[0].isdigit():
                continue
            if any(ch in tok for ch in "()+²√θ=,:{}^\\") and not re.search(r"\.lean:\d+$", tok):
                continue  # a formula, or a toolchain name such as `leanprover/lean4:v4.33.0`
            if os.path.basename(tok.rstrip("/")) in external:
                continue
            if not ("/" in tok or tok.endswith(exts)):
                continue
            if tok.count("/") == 1 and tok.endswith("/") is False and not tok.endswith(exts) and "." not in tok.split("/")[-1] and tok.split("/")[0] in ("n_i", "c_i", "1", "θ", "d", "n", "c"):
                continue  # a fraction, not a path
            if re.fullmatch(r"[A-Za-z0-9_.]+/[A-Za-z0-9_.]+", tok) and tok.split("/")[0] not in ("docs", "notes", "scripts", "lean", "StackedSVD", "paper", "Prob", "RMT", "LinAlg", "RankR", "SVDStack", "StackSVD", "General", "Vendor", "MLEMarginal", "Het", "SingleWeight", "archive", "agent_reports", "audit_packets", "discoveries", "numeric", "heteroedge_elementary", ".github", "workflows", "COLT83", "src", "StatsMLlib"):
                continue  # `a/b` outside any known folder: a ratio or a quotient, not a path
            where = resolve(tok, docdir, all_files, tree_set, on_disk=True)
            if where is None and withheld:
                where = resolve(tok, docdir, withheld, withheld_set)
                if where is not None:
                    n_withheld += 1
                    if mode == "--list":
                        print(f"{doc}:{ln}: `{tok}` -> {where} (withheld from this release)")
                    continue
            if where is None:
                around = prev[-120:] + " " + line[max(0, m.start() - 120):m.end() + 60]
                if skip_words.search(around):
                    continue
                n_bad += 1; bad.append(f"{doc}:{ln}: `{tok}`")
            else:
                n_ok += 1
                if mode == "--list":
                    print(f"{doc}:{ln}: `{tok}` -> {where}")
        prev = line

for b in bad:
    print("check_paths: missing", b)
if n_bad:
    print(f"check_paths: FAIL - {n_bad} path(s) do not resolve, {n_ok} do.")
    sys.exit(1)
tail = f", {n_withheld} name files withheld from this release (scripts/release_withheld.txt)" if withheld else ""
print(f"check_paths: OK - {n_ok} paths resolve{tail}, 0 missing.")
