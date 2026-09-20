#!/bin/bash
# check_layering.sh - declaration-level check that the Layer 2 Gaussian discharges do not
# use a Layer 1 theorem that consumes the hypothesis structure they exist to prove.
#
# `RMT/Het/MPhet.lean`, `RMT/Het/Split.lean` and `RMT/Het/Simplicity.lean` import
# `StackSVDWeighted.lean` at the file level, because they need the type `HeteroLaw` to
# state what they prove. This script runs `scripts/LayerAudit.lean`, which walks the proof
# term of each of the three Gaussian discharges (`singleTableLaw_of_gaussian`,
# `heteroLaw_of_gaussian`, `heteroEdge_of_gaussian`) and confirms that none of them calls a
# Layer 1 theorem such as `thm_stacksvd_weighted` or `thm_theta_est`. Two criteria: a fixed
# list of 13 Layer 1 theorem names (`FORBIDDEN` lines), and a list-free one, any constant in
# the closure whose type takes a `SingleTableLaw` or a `HeteroLaw` as a hypothesis
# (`CONSUMER` lines; the structures' own projections are exempt). The driver self-tests the
# second criterion on five known consumers before the audit. File-level imports and
# declaration-level use are different questions; this answers the second one.
#
# Usage (from anywhere):
#   scripts/check_layering.sh            gate: exit 0 iff 3 endpoints, 0 forbidden hits
#   scripts/check_layering.sh --raw      print the driver output only, no verdict
#   scripts/check_layering.sh --help
#
# The driver reads the imported environment only, so it reports the state of the oleans in
# `lean/.lake/build`, not the state of the sources (same caveat as check_axioms.sh;
# this script does not check freshness).
#
# Exit codes: 0 pass, 1 a forbidden constant turned up or the endpoint count is wrong,
# 2 usage, environment, or Lean error (including a kill by the 1500 s timeout).
set -eo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$here")"
proj="$root/lean"
tools="${LEAN_TOOLS:-}"   # helper scripts of the development server; absent elsewhere
env_sh="$tools/env.sh"
driver="$here/LayerAudit.lean"
mode="${1:-gate}"

case "$mode" in
  gate|--raw) ;;
  --help|-h) sed -n '2,19p' "$here/check_layering.sh"; exit 0 ;;
  *) echo "check_layering.sh: unknown option '$mode'" >&2; exit 2 ;;
esac

[ -d "$proj" ] || { echo "check_layering.sh: no lake project at $proj" >&2; exit 2; }
[ -f "$driver" ] || { echo "check_layering.sh: no $driver" >&2; exit 2; }
# The server keeps its toolchain and cache off the home filesystem and selects them in
# env.sh; on any other machine the elan-managed `lake` on PATH is used.
if [ -f "$env_sh" ]; then
  # shellcheck disable=SC1090
  source "$env_sh"
elif ! command -v lake >/dev/null; then
  echo "check_layering.sh: lake is not on PATH (install elan, https://github.com/leanprover/elan)" >&2
  exit 2
fi

# ---------------------------------------------------------------- run the driver
out="$(mktemp)"
trap 'rm -f "$out"' EXIT

lean_status=0
# GNU `timeout` caps the run at 1500 s where it exists (Linux); macOS has none by default.
if command -v timeout >/dev/null; then tmo="timeout 1500"; else tmo=""; fi
( cd "$proj" && $tmo lake env lean -j 2 "$driver" ) >"$out" 2>&1 || lean_status=$?

cat "$out"

if [ "$lean_status" -eq 124 ]; then
  echo "check_layering: FAIL - lake env lean did not finish inside 1500 s." >&2
  exit 2
elif [ "$lean_status" -ne 0 ]; then
  echo "check_layering: FAIL - lake env lean exited $lean_status (a Lean error, or an" \
    "unresolved endpoint or forbidden name; see LayerAudit.lean's throwError)." >&2
  exit 2
fi

if [ "$mode" = "--raw" ]; then
  exit 0
fi

# ---------------------------------------------------------------- verdict
summary="$(grep '^LAYERAUDIT_END' "$out" || true)"
if [ -z "$summary" ]; then
  echo "check_layering: FAIL - no LAYERAUDIT_END line in the driver output." >&2
  exit 2
fi

read -r _ n_endpoints n_forbidden <<< "$summary"

if [ "$n_endpoints" = "3" ] && [ "$n_forbidden" = "0" ]; then
  echo "check_layering: OK - 3 endpoints audited, 0 forbidden hits."
  exit 0
else
  echo "check_layering: FAIL - $n_endpoints endpoint(s) audited, $n_forbidden forbidden" \
    "hit(s) (expected 3 and 0)." >&2
  exit 1
fi
