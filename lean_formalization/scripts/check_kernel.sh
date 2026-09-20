#!/bin/bash
# check_kernel.sh - gate 7: replay every module of the package through the Lean kernel a
# second time, with the toolchain's own `leanchecker` (Lean v4.28 and later; before that the
# separate lean4checker, now archived).
#
# What it checks. `lake build` elaborates a file and adds each declaration to the environment
# through the kernel once. `leanchecker M` reads the compiled `.olean` of `M`, imports the
# oleans of the imports of `M`, and adds every constant of `M` to that environment through
# the kernel again (`Lean.Environment.replay`), with no elaborator, no tactic and no `unsafe`
# code in between. It catches environment hacking (a declaration that reached the olean
# without a kernel check, for example through `addDeclWithoutChecking` or a `set_option
# debug.skipKernelTC`) and a corrupt olean. It trusts the imports as loaded: the gate replays
# every module of ours, so each module of ours is re-checked when the gate reaches it, and
# Mathlib and StatsMLlib are trusted as compiled. It is not an external verifier: the kernel
# that re-checks is the one of the same toolchain (lean4lean would be the independent one).
#
# How it runs on the development server. Import time dominates: each replay imports the
# Mathlib closure once, about 5000 oleans, 7 s on the root-disk copy that `lean-local.sh sync`
# keeps at $LEAN_LOCAL_LIBS (default /tmp/$USER-lean-libs) and 5 to 10 min over NFS. So the gate runs leanchecker directly on the
# copy (`lean-local.sh path`) after it has checked that every `.olean` of ours in the copy is
# byte-identical (`cmp`) to the build output in `lean/.lake/build/lib/lean`, so what it
# replays is what `lake build` wrote; a difference means "run lean-local.sh sync" and exit 2.
# On a machine without the copy it falls back to the search path of `lake env` (correct,
# slow over NFS, fine on a local disk).
#
# Concurrency. Both modes replay KERNEL_CHECK_JOBS modules at once (default 3, about 17 GB;
# each replay imports the Mathlib closure and holds about 5 GB). leanchecker has no `-j`: it
# starts one task per matched module and the Lean runtime sizes the pool by the CPU count, so
# the 182 tasks of `leanchecker StackedSVD` reached 231 GB (2026-09-07). Mode 1, with the
# nprocs shim of lake-build-capped.sh, caps that pool. Mode 2, the portable one, runs
# scripts/KernelReplay.lean instead: a copy of leanchecker whose own pool is bounded, so the
# cap needs no shim and works on any machine.
#
# Usage (from anywhere):
#   scripts/check_kernel.sh                                    gate: every module of the package
#   scripts/check_kernel.sh StackedSVD.Defs StackedSVD.RankR   only these module-name prefixes
#   KERNEL_CHECK_JOBS=2 scripts/check_kernel.sh                a smaller pool
#   scripts/check_kernel.sh --help
#
# Output: one `replaying <module>` line per module (both modes), then the verdict
#   check_kernel: OK - <n> modules replayed through the kernel, 0 problems (<s> s, <j> workers)
# and, for the default target, a comparison of <n> with the number of oleans of the package.
#
# Exit codes: 0 pass, 1 the replay reported a problem or the module count is off, 2 usage or
# environment error (no toolchain, no oleans, a stale copy).
set -eo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$here")"
proj="$root/lean"
# LEAN_TOOLS: the helper scripts of the development server (env.sh, lean-local.sh, the nprocs
# shim); unset or absent elsewhere, and then the script uses `lake` from PATH with no cap.
tools="${LEAN_TOOLS:-}"
env_sh="$tools/env.sh"
local_sh="$tools/lean-local.sh"
shim="$tools/shim/nprocs_shim.so"
libs="${LEAN_LOCAL_LIBS:-/tmp/$USER-lean-libs}"   # the root-disk copy of lean-local.sh
jobs="${KERNEL_CHECK_JOBS:-3}"

case "${1:-}" in
  --help|-h) sed -n '2,45p' "$here/check_kernel.sh"; exit 0 ;;
  -*) echo "check_kernel.sh: unknown option '$1'" >&2; exit 2 ;;
esac
targets=("$@")
[ ${#targets[@]} -gt 0 ] || targets=(StackedSVD)

[ -d "$proj/.lake/build/lib/lean/StackedSVD" ] || {
  echo "check_kernel.sh: no build output at $proj/.lake/build/lib/lean (run lake build first)" >&2
  exit 2
}
if [ -f "$env_sh" ]; then
  # shellcheck disable=SC1090
  source "$env_sh"
elif ! command -v lake >/dev/null; then
  echo "check_kernel.sh: lake is not on PATH (install elan, https://github.com/leanprover/elan)" >&2
  exit 2
fi

# ---------------------------------------------------------------- the toolchain
sysroot="$(cd "$proj" && lean --print-prefix 2>/dev/null)" || sysroot=""
[ -n "$sysroot" ] && [ -x "$sysroot/bin/leanchecker" ] || {
  echo "check_kernel.sh: no leanchecker in the toolchain of $proj (needs Lean v4.28 or later)" >&2
  exit 2
}
export LEAN_SYSROOT="$sysroot"   # leanchecker asks `lean --print-prefix` otherwise, which needs a toolchain in the cwd
checker="$sysroot/bin/leanchecker"

# ---------------------------------------------------------------- the search path
build_lib="$proj/.lake/build/lib/lean"
copy=""
if [ -x "$local_sh" ] && [ -f "$libs/StackedSVD/.complete" ]; then
  copy="$libs/StackedSVD/lib/lean"
  n_cmp=0; n_bad=0
  while IFS= read -r -d '' f; do
    rel="${f#"$build_lib"/}"
    n_cmp=$((n_cmp + 1))
    if ! cmp -s "$f" "$copy/$rel"; then
      n_bad=$((n_bad + 1))
      [ $n_bad -le 5 ] && echo "check_kernel.sh: copy differs or is missing: $rel" >&2
    fi
  done < <(find "$build_lib/StackedSVD" "$build_lib/StackedSVD.olean" -name '*.olean' -print0)
  if [ $n_bad -gt 0 ]; then
    echo "check_kernel.sh: $n_bad of $n_cmp oleans of the package differ between the build and the" \
         "root-disk copy; run '$local_sh sync' after the build, then rerun." >&2
    exit 2
  fi
  export LEAN_PATH="$("$local_sh" path)"
  where="the root-disk copy ($n_cmp oleans of the package byte-identical to the build)"
else
  export LEAN_PATH="$(cd "$proj" && lake env printenv LEAN_PATH)"
  where="the lake search path"
fi

# ---------------------------------------------------------------- run
if [ -f "$shim" ] && command -v taskset >/dev/null; then
  # Mode 1: leanchecker, with the nprocs shim of lake-build-capped.sh around its task pool.
  engine="leanchecker of $(basename "$sysroot")"
  run=(env "FAKE_NPROCS=$jobs" "LD_PRELOAD=$shim" taskset -c "0-$((jobs - 1))" \
       "$checker" -v "${targets[@]}")
else
  # Mode 2, portable: KernelReplay.lean, a copy of leanchecker whose pool is bounded by -j.
  # The first -j sizes the thread pool of the Lean runtime. It matters: at the CPU count (64
  # on the development server) the reserved address space passes 40 GB and a `ulimit -v
  # 40000000` kills the run. The second -j sizes the pool of replays.
  [ -f "$here/KernelReplay.lean" ] || {
    echo "check_kernel.sh: no $here/KernelReplay.lean (needed without the nprocs shim)" >&2
    exit 2
  }
  engine="KernelReplay.lean on $(basename "$sysroot") (portable mode, no nprocs shim)"
  run=("$sysroot/bin/lean" -j "$jobs" --run "$here/KernelReplay.lean" -v -j "$jobs" "${targets[@]}")
  echo "check_kernel: no nprocs shim; portable mode, $jobs replays at once (about 5 GB each)" >&2
fi
out="$(mktemp)"
trap 'rm -f "$out"' EXIT
start=$(date +%s)
echo "check_kernel: $engine, $jobs workers, oleans from $where"
echo "check_kernel: targets ${targets[*]}"
status=0
"${run[@]}" 2>&1 | tee "$out" || status=${PIPESTATUS[0]}
secs=$(( $(date +%s) - start ))
n_rep=$(grep -c '^replaying ' "$out" || true)

if [ "$status" != "0" ] || grep -q 'leanchecker found a problem' "$out"; then
  echo "check_kernel: FAIL - the replay exited $status after $n_rep replays ($secs s); see the lines above."
  exit 1
fi
if [ "${targets[*]}" = "StackedSVD" ]; then
  n_ole=$(find "$build_lib/StackedSVD" "$build_lib/StackedSVD.olean" -name '*.olean' | wc -l)
  if [ "$n_rep" != "$n_ole" ]; then
    echo "check_kernel: FAIL - $n_rep modules replayed but the package has $n_ole oleans ($secs s)."
    exit 1
  fi
  echo "check_kernel: OK - $n_rep modules replayed through the kernel, 0 problems; every one of the" \
       "$n_ole oleans of the package ($secs s, $jobs workers)."
else
  [ "$n_rep" -gt 0 ] || { echo "check_kernel: FAIL - nothing replayed ($secs s)."; exit 1; }
  echo "check_kernel: OK - $n_rep modules replayed through the kernel, 0 problems ($secs s, $jobs workers)."
fi
