#!/bin/bash
# check_sorries.sh - compare SORRIES.md with the `sorry` sites in the tree.
#
# Usage (from anywhere):
#   scripts/check_sorries.sh            gate: exit 0 if the table matches the tree, 1 if not
#   scripts/check_sorries.sh --sites    print `file:line<TAB>declaration` for every site
#   scripts/check_sorries.sh --tracked  print the declaration of every SORRIES.md row
#   scripts/check_sorries.sh --help
#
# Scope: lean/StackedSVD/**/*.lean, without Vendor/ (third party code).
# A site is a `sorry` term or tactic. Line comments, block comments (nested), doc
# comments, string literals and character literals are blanked before the search, so
# the word `sorry` inside a comment or a message is not a site.
#
# Exit codes: 0 match (or listing mode), 1 mismatch, 2 usage or environment error.
set -eo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$here")"
mode="${1:-gate}"

case "$mode" in
  gate|--sites|--tracked) ;;
  --help|-h) sed -n '2,15p' "$here/check_sorries.sh"; exit 0 ;;
  *) echo "check_sorries.sh: unknown option '$mode'" >&2; exit 2 ;;
esac

[ -f "$root/docs/SORRIES.md" ] || { echo "check_sorries.sh: no docs/SORRIES.md under $root" >&2; exit 2; }
[ -d "$root/lean/StackedSVD" ] || { echo "check_sorries.sh: no package dir" >&2; exit 2; }

python3 - "$root" "$mode" <<'PYEOF'
import os, re, sys

root, mode = sys.argv[1], sys.argv[2]
pkg = os.path.join(root, "lean", "StackedSVD")
md_path = os.path.join(root, "docs", "SORRIES.md")

# ---------------------------------------------------------------- comment stripper
CHAR_RE = re.compile(r"'(\\[^']{0,8}|[^'\\])'")

def blank(text):
    """Replace comment, string and character literal content with spaces.
    Newlines survive, so line numbers do not move."""
    out, i, n, depth = [], 0, len(text), 0
    while i < n:
        ch = text[i]
        nxt = text[i + 1] if i + 1 < n else ""
        if depth > 0:                                   # inside /- ... -/
            if ch == "/" and nxt == "-":
                depth += 1; out.append("  "); i += 2; continue
            if ch == "-" and nxt == "/":
                depth -= 1; out.append("  "); i += 2; continue
            out.append("\n" if ch == "\n" else " "); i += 1; continue
        if ch == "/" and nxt == "-":
            depth = 1; out.append("  "); i += 2; continue
        if ch == "-" and nxt == "-":                    # -- to end of line
            j = text.find("\n", i)
            j = n if j < 0 else j
            out.append(" " * (j - i)); i = j; continue
        if ch == '"':                                   # string literal
            out.append(" "); i += 1
            while i < n:
                c = text[i]
                if c == "\\" and i + 1 < n:
                    out.append("  "); i += 2; continue
                if c == '"':
                    out.append(" "); i += 1; break
                out.append("\n" if c == "\n" else " "); i += 1
            continue
        if ch == "'" and (i == 0 or not (text[i - 1].isalnum() or text[i - 1] in "_'")):
            m = CHAR_RE.match(text, i)                  # character literal, not a prime
            if m:
                out.append(" " * (m.end() - m.start())); i = m.end(); continue
        out.append(ch); i += 1
    return "".join(out)

# ---------------------------------------------------------------- site detection
SORRY_RE = re.compile(r"(?<![\w'.₀-₟])sorry(?![\w'!?])")
DECL_RE = re.compile(
    r"^\s*(?:@\[[^\]]*\]\s*)*"
    r"(?:(?:private|protected|public|noncomputable|partial|unsafe|scoped|local|nonrec|meta)\s+)*"
    r"(theorem|lemma|def|instance|abbrev|example|structure|class|inductive)\b\s*(.*)$")
NAME_RE = re.compile(r"[^\s:({\[⦃⟨⟪]+")

def decl_of(lines, idx):
    """Nearest preceding declaration header at or above line index `idx`."""
    for k in range(idx, -1, -1):
        m = DECL_RE.match(lines[k])
        if not m:
            continue
        kw, rest = m.group(1), m.group(2).strip()
        nm = NAME_RE.match(rest)
        name = nm.group(0).rstrip(",") if nm else ""
        if not name or name in ("where", "extends"):
            name = "<anonymous %s line %d>" % (kw, k + 1)
        return name
    return "<no enclosing declaration>"

lean_files = []
for dirpath, dirnames, filenames in os.walk(pkg):
    dirnames[:] = [d for d in dirnames if d != "Vendor"]
    for f in sorted(filenames):
        if f.endswith(".lean"):
            lean_files.append(os.path.join(dirpath, f))
lean_files.sort()

sites = []          # (relpath, line, declaration)
for path in lean_files:
    rel = os.path.relpath(path, pkg)
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    stripped = blank(text)
    lines = stripped.split("\n")
    for m in SORRY_RE.finditer(stripped):
        ln = stripped.count("\n", 0, m.start())
        sites.append((rel, ln + 1, decl_of(lines, ln)))

# ---------------------------------------------------------------- SORRIES.md rows
def canon_file(cell):
    cand = cell
    for _ in range(2):
        if os.path.isfile(os.path.join(pkg, cand)):
            return cand
        if cand.startswith("StackedSVD/"):
            cand = cand[len("StackedSVD/"):]
        else:
            break
    return cand

rows = []           # (canonical relpath, declaration as written)
with open(md_path, encoding="utf-8") as fh:
    for raw in fh:
        s = raw.strip()
        if not s.startswith("|"):
            continue
        cells = [c.strip() for c in s.strip("|").split("|")]
        if len(cells) < 2:
            continue
        if set(cells[0]) <= set("-: "):
            continue
        f = cells[0].strip("`").strip()
        d = cells[1].strip("`").strip()
        if f.lower() == "file" and d.lower() == "declaration":
            continue
        rows.append((canon_file(f), d))

# ---------------------------------------------------------------- listing modes
if mode == "--sites":
    for rel, ln, decl in sites:
        print("lean/StackedSVD/%s:%d\t%s" % (rel, ln, decl))
    sys.exit(0)
if mode == "--tracked":
    for _, d in rows:
        print(d)
    sys.exit(0)

# ---------------------------------------------------------------- compare
def key(f, d):
    return (f, d.split(".")[-1])

site_keys = {}
for rel, ln, decl in sites:
    site_keys.setdefault(key(rel, decl), []).append((rel, ln, decl))
row_keys = {}
for f, d in rows:
    row_keys.setdefault(key(f, d), []).append((f, d))

missing_site = [k for k in row_keys if k not in site_keys]     # a row with no `sorry`
missing_row = [k for k in site_keys if k not in row_keys]      # a `sorry` with no row

if not missing_site and not missing_row:
    print("check_sorries: OK - %d sorry sites in %d declarations, %d rows, all matched."
          % (len(sites), len(site_keys), len(rows)))
    sys.exit(0)

print("check_sorries: FAIL - SORRIES.md does not match the tree.")
if missing_row:
    print("\n  Untracked sorry sites (add a row to SORRIES.md):")
    for k in sorted(missing_row):
        for rel, ln, decl in site_keys[k]:
            print("    lean/StackedSVD/%s:%d  %s" % (rel, ln, decl))
if missing_site:
    print("\n  Rows with no sorry in the tree (delete the row, or the file moved):")
    for k in sorted(missing_site):
        for f, d in row_keys[k]:
            print("    StackedSVD/%s  %s" % (f, d))
print("\n  Sites found: %d. Rows read: %d." % (len(sites), len(rows)))
sys.exit(1)
PYEOF
