#!/bin/bash
# make_audit_pack_header.sh - generate the dated build-and-gates header (once for notes/archive/AUDIT_PACK_V2.md, now docs/AUDIT_DOC.md Appendix A).
#
# One script produces every dated or state-dependent line of the audit pack header,
# so no date, commit or build result in the pack is hand written.
#
# It prints, in order:
#   1. the UTC time it started (date -u)
#   2. the git commit (rev-parse HEAD, log --oneline -1) and the worktree state
#   3. the Lean toolchain (lean/lean-toolchain)
#   4. a full build: $LEAN_TOOLS/lake-build-capped.sh 12 (the development server), else lake build
#      (last line of the build log, the exit code, and how many of the jobs were
#      recompiled as opposed to replayed from the existing .lake cache)
#   5. scripts/check_sorries.sh          (full output and exit code)
#   6. scripts/check_axioms.sh           (full output and exit code)
#   7. scripts/check_root_imports.py     (full output and exit code)
#   8. scripts/check_theorems_sigs.py    (full output and exit code)
#   9. the UTC time it finished
#
# Usage (from anywhere):
#   scripts/make_audit_pack_header.sh              build and print the header
#   scripts/make_audit_pack_header.sh --no-build   skip step 4 (for a dry run only)
#   AUDIT_PACK_JOBS=6 scripts/make_audit_pack_header.sh    cap the build at 6 cores
#
# AUDIT_PACK_JOBS sets the core count of the capped build (default 12). Use it when the
# session has a lower thread budget; the printed header names the count it used.
#
# The script always exits 0. Read the recorded exit codes: a non-zero code there
# means the pack must not claim that the library builds.
set -eo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$here")"
proj="$root/lean"
tools="${LEAN_TOOLS:-}"   # helper scripts of the development server; absent elsewhere
env_sh="$tools/env.sh"
capped="$tools/lake-build-capped.sh"
jobs="${AUDIT_PACK_JOBS:-12}"
do_build=1
[ "${1:-}" = "--no-build" ] && do_build=0

# shellcheck disable=SC1090
source "$env_sh" >/dev/null 2>&1 || true

echo "generated-by : scripts/make_audit_pack_header.sh"
echo "start (UTC)  : $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
echo "host         : $(hostname)"
echo ""
echo "--- git ---"
echo "HEAD         : $(cd "$root" && git rev-parse HEAD)"
echo "log -1       : $(cd "$root" && git log --oneline -1)"
dirty="$(cd "$root" && git status --porcelain | wc -l)"
echo "worktree     : $dirty file(s) modified or untracked at header time"
(cd "$root" && git status --porcelain | sed 's/^/               /' | head -20)
echo ""
echo "--- toolchain ---"
echo "lean-toolchain: $(cat "$proj/lean-toolchain")"
echo "lake version  : $( (cd "$proj" && lake --version) 2>&1 | head -1)"
echo ""
echo "--- build: lake-build-capped.sh $jobs ---"
if [ "$do_build" = "1" ]; then
  log="$(mktemp)"
  rc=0
  ( cd "$proj" && "$capped" "$jobs" ) >"$log" 2>&1 || rc=$?
  echo "last line    : $(tail -n 1 "$log")"
  echo "exit code    : $rc"
  echo "log lines    : $(wc -l <"$log")"
  echo "errors       : $(grep -c '^error' "$log" || true)"
  echo "jobs rebuilt : $(grep -c '\] Built ' "$log" || true)"
  echo "jobs replayed: $(grep -c 'Replayed' "$log" || true)  (cached traces re-checked, not recompiled)"
  rm -f "$log"
else
  echo "(skipped: --no-build)"
fi
echo ""
echo "--- gate: scripts/check_sorries.sh ---"
rc=0
bash "$here/check_sorries.sh" || rc=$?
echo "exit code    : $rc"
echo ""
echo "--- gate: scripts/check_axioms.sh ---"
rc=0
bash "$here/check_axioms.sh" || rc=$?
echo "exit code    : $rc"
echo ""
echo "--- gate: scripts/check_root_imports.py ---"
rc=0
python3 "$here/check_root_imports.py" || rc=$?
echo "exit code    : $rc"
echo ""
echo "--- gate: scripts/check_theorems_sigs.py docs/THEOREMS.md ---"
rc=0
python3 "$here/check_theorems_sigs.py" "$root/docs/THEOREMS.md" || rc=$?
echo "exit code    : $rc"
echo ""
echo "end (UTC)    : $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
