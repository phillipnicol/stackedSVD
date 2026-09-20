#!/bin/bash
# check_axioms.sh - axiom gate for the StackedSVD library.
#
# It runs `scripts/AxiomAudit.lean` against the built oleans and fails when a
# declaration of the `StackedSVD` namespace depends on anything but
# `propext`, `Classical.choice` and `Quot.sound`.
#
# `sorryAx` is the one conditional case. It passes when the declaration has a row in
# SORRIES.md, and when a declaration only inherits it from such a row (the driver
# reports the source of every `sorry` on a `SORRYDEPS` line). Any other `sorryAx`
# fails.
#
# Usage (from anywhere):
#   scripts/check_axioms.sh            gate
#   scripts/check_axioms.sh --raw      print the driver output, no gate
#   scripts/check_axioms.sh --verbose  gate, and print every declaration
#   scripts/check_axioms.sh --help
#
# Environment:
#   CHECK_AXIOMS_BUILD=1   build first (`lake-build-capped.sh 12` on the server, else `lake build`). Default is no build:
#                          the script only warns when the oleans look stale.
#   AXIOM_AUDIT_VENDOR=1   also audit `StackedSVD/Vendor/` (third party backports).
#
# CAUTION: without CHECK_AXIOMS_BUILD=1 the gate reads the oleans, not the sources.
# A declaration edited since the last build is audited in its old form, and a new
# declaration is absent. Read the freshness warning that the script prints.
#
# Exit codes: 0 pass, 1 axiom violation, 2 usage or environment or build error.
set -eo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$here")"
proj="$root/lean"
pkg="$proj/StackedSVD"
# LEAN_TOOLS: the helper scripts of the development server (env.sh, lake-build-capped.sh);
# unset or absent elsewhere, and then the script uses `lake` from PATH.
tools="${LEAN_TOOLS:-}"
env_sh="$tools/env.sh"
capped="$tools/lake-build-capped.sh"
mode="${1:-gate}"

case "$mode" in
  gate|--raw|--verbose) ;;
  --help|-h) sed -n '2,28p' "$here/check_axioms.sh"; exit 0 ;;
  *) echo "check_axioms.sh: unknown option '$mode'" >&2; exit 2 ;;
esac

[ -d "$pkg" ] || { echo "check_axioms.sh: no package dir at $pkg" >&2; exit 2; }
# The server keeps its toolchain and cache off the home filesystem and selects them in
# env.sh; on any other machine the elan-managed `lake` on PATH is used.
if [ -f "$env_sh" ]; then
  # shellcheck disable=SC1090
  source "$env_sh"
elif ! command -v lake >/dev/null; then
  echo "check_axioms.sh: lake is not on PATH (install elan, https://github.com/leanprover/elan)" >&2
  exit 2
fi

# ---------------------------------------------------------------- freshness report
# `python3` does the timestamps (GNU `find -printf` and `date -d` do not exist on macOS).
extremum() {  # extremum <dir> <suffix> <max|min>: "<epoch> <path>" of the newest or oldest file
  python3 - "$1" "$2" "$3" <<'PY'
import os, sys
d, suf, which = sys.argv[1:4]
best = None
for dp, _, fs in os.walk(d):
    for f in fs:
        if f.endswith(suf):
            p = os.path.join(dp, f); t = os.path.getmtime(p)
            if best is None or (t > best[0] if which == "max" else t < best[0]): best = (t, p)
print(f"{best[0]} {best[1]}" if best else "")
PY
}
newest_src="$(extremum "$pkg" .lean max)"
oldest_olean="$(extremum "$proj/.lake/build/lib/lean" .olean min 2>/dev/null || true)"
fmt() {
  [ -n "$1" ] || { echo "none"; return; }
  python3 -c 'import sys, time; print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(float(sys.argv[1]))))' "${1%% *}"
}

if [ "${CHECK_AXIOMS_BUILD:-0}" = "1" ]; then
  if pgrep -x lake >/dev/null; then
    echo "check_axioms: another lake process runs; refusing to build." >&2
    exit 2
  fi
  if [ -x "$capped" ]; then
    echo "check_axioms: CHECK_AXIOMS_BUILD=1, running $capped 12"
    ( cd "$proj" && "$capped" 12 ) || { echo "check_axioms: build failed." >&2; exit 2; }
  else
    echo "check_axioms: CHECK_AXIOMS_BUILD=1, running lake build"
    ( cd "$proj" && lake build ) || { echo "check_axioms: build failed." >&2; exit 2; }
  fi
else
  echo "check_axioms: WARNING - no build. The audit reads the oleans in .lake/build."
  echo "  newest source: $(fmt "$newest_src")  ${newest_src#* }"
  echo "  oldest olean : $(fmt "$oldest_olean")  ${oldest_olean#* }"
  # per file: a source newer than its own olean is stale (audit_mechanical 2026-08-31, finding 6:
  # the old check compared the newest source with the oldest olean and always fired)
  miss=0; stale=0
  while IFS= read -r f; do
    rel="${f#"$pkg"/}"; rel="${rel%.lean}"
    ol="$proj/.lake/build/lib/lean/StackedSVD/$rel.olean"
    if [ ! -f "$ol" ]; then
      echo "  NO OLEAN: StackedSVD/$rel.lean"; miss=$((miss + 1))
    elif [ "$f" -nt "$ol" ]; then
      echo "  STALE: StackedSVD/$rel.lean is newer than its olean"; stale=$((stale + 1))
    fi
  done < <(find "$pkg" -name '*.lean' | sort)
  [ "$miss" -gt 0 ] && echo "  $miss source file(s) are not in the build; their declarations are not audited."
  [ "$stale" -gt 0 ] && echo "  $stale source file(s) are newer than their olean. Set CHECK_AXIOMS_BUILD=1 for a real gate."
fi

# ---------------------------------------------------------------- run the driver
out="$(mktemp)"; err="$(mktemp)"
trap 'rm -f "$out" "$err"' EXIT
rc=0
( cd "$proj" && lake env lean -j 3 "$here/AxiomAudit.lean" ) >"$out" 2>"$err" || rc=$?
if [ "$rc" -ne 0 ] && pgrep -x lake >/dev/null; then
  # A concurrent `lake build` rewrites the oleans, so a read can fail. Try once more.
  echo "check_axioms: the driver failed while a lake build runs; one retry."
  rc=0
  ( cd "$proj" && lake env lean -j 3 "$here/AxiomAudit.lean" ) >"$out" 2>"$err" || rc=$?
fi
if [ "$rc" -ne 0 ]; then
  echo "check_axioms: the driver failed (exit $rc)." >&2
  sed -n '1,40p' "$err" >&2
  sed -n '1,20p' "$out" >&2
  pgrep -x lake >/dev/null && echo "check_axioms: a lake build runs; the oleans move under the driver." >&2
  exit 2
fi
[ -s "$err" ] && { echo "check_axioms: driver messages:"; sed -n '1,20p' "$err"; }

if [ "$mode" = "--raw" ]; then cat "$out"; exit 0; fi

# ---------------------------------------------------------------- tracked sorries
tracked="$(mktemp)"; trap 'rm -f "$out" "$err" "$tracked"' EXIT
if ! bash "$here/check_sorries.sh" --tracked >"$tracked"; then
  echo "check_axioms: could not read SORRIES.md." >&2; exit 2
fi
if ! bash "$here/check_sorries.sh" >/dev/null 2>&1; then
  echo "check_axioms: NOTE - check_sorries.sh does not pass, so the tracked set may be wrong."
fi

# ---------------------------------------------------------------- verdict
python3 - "$out" "$tracked" "$mode" <<'PYEOF'
import re, sys

out_path, tracked_path, mode = sys.argv[1], sys.argv[2], sys.argv[3]
ALLOWED = {"propext", "Classical.choice", "Quot.sound"}

tracked = [l.strip() for l in open(tracked_path, encoding="utf-8") if l.strip()]

AUX = re.compile(r"\.(?:proof|match|eq|omega)_\d+$|\._@.*$|\.\d+$")
PRIV = re.compile(r"^_private\.[^.]+(?:\.[^.]+)*?\.\d+\.")

def normalize(name):
    n = PRIV.sub("", name)
    prev = None
    while prev != n:
        prev, n = n, AUX.sub("", n)
    return n

def is_tracked(name):
    n = normalize(name)
    return any(n == t or n.endswith("." + t) for t in tracked)

PRIVATE_MARK = " (private)"

def split_record(payload):
    name, _, rest = payload.partition("|")
    name = name.strip()
    priv = name.endswith(PRIVATE_MARK)
    if priv:
        name = name[: -len(PRIVATE_MARK)].strip()
    items = [x for x in rest.strip().split(",") if x and x != "-"]
    return name, priv, items

decls, deps, module_of, cur_mod = [], {}, {}, "?"
end = None
for line in open(out_path, encoding="utf-8"):
    line = line.rstrip("\n")
    if line.startswith("MODULE "):
        cur_mod = line[7:].strip()
    elif line.startswith("DECL "):
        name, priv, axs = split_record(line[5:])
        decls.append((name, axs))
        module_of[name] = cur_mod + (" private" if priv else "")
    elif line.startswith("SORRYDEPS "):
        name, _priv, src = split_record(line[10:])
        deps[name] = src
    elif line.startswith("AXIOMAUDIT_END"):
        end = line

if end is None:
    print("check_axioms: FAIL - the driver output is truncated (no AXIOMAUDIT_END).")
    sys.exit(1)

bad_axiom, bad_sorry, ok_sorry, inherited = [], [], [], []
for name, axs in decls:
    extra = [a for a in axs if a not in ALLOWED and a != "sorryAx"]
    if extra:
        bad_axiom.append((name, extra))
    if "sorryAx" in axs:
        srcs = deps.get(name, [])
        untracked = [s for s in srcs if not is_tracked(s)]
        if not srcs:
            bad_sorry.append((name, ["<no source found>"]))
        elif untracked:
            bad_sorry.append((name, untracked))
        elif is_tracked(name):
            ok_sorry.append(name)
        else:
            inherited.append((name, srcs))

if mode == "--verbose":
    for name, axs in decls:
        print("  %-70s %s" % (name, ",".join(axs) if axs else "-"))

print("check_axioms: %s declarations audited (%s)." % (len(decls), end.split(" ", 1)[1]))
if ok_sorry:
    print("  tracked sorry (SORRIES.md row):")
    for n in sorted(ok_sorry):
        print("    %s" % n)
if inherited:
    print("  inherits a tracked sorry (no row needed):")
    for n, s in sorted(inherited):
        print("    %s  <- %s" % (n, ",".join(s)))

if not bad_axiom and not bad_sorry:
    print("check_axioms: OK - no axiom outside {propext, Classical.choice, Quot.sound}"
          " and no untracked sorry.")
    sys.exit(0)

print("check_axioms: FAIL")
for n, a in sorted(bad_axiom):
    print("  AXIOM  %s  depends on %s  [%s]" % (n, ",".join(a), module_of.get(n, "?")))
for n, s in sorted(bad_sorry):
    print("  SORRY  %s  untracked sorry from %s  [%s]" % (n, ",".join(s), module_of.get(n, "?")))
sys.exit(1)
PYEOF
