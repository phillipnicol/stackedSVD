/-
Copyright (c) 2023 Kim Morrison. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Kim Morrison, Sebastian Ullrich

This file is a copy of `src/lean/LeanChecker.lean` of the Lean 4 toolchain (v4.33.0), the
source of the `leanchecker` executable, with the changes listed below. The original license
and the original authors stay as they are.

Why the copy. `leanchecker` starts one `IO.asTask` per matched module at once. The Lean task
pool then runs as many as the machine has cores, and each running replay holds about 5 GB
(it imports the Mathlib closure). On a 16 core machine `leanchecker StackedSVD` peaks near
80 GB. `leanchecker` has no `-j`, and its targets are module name prefixes, so a shell loop
over one module at a time does not work either: three modules of this package have both an
olean and a folder of children.

What changed against the original:
1. A bounded worker pool. At most `-j N` replays are alive at once (default 3, about 15 GB).
   The pool is a sliding window: it starts N tasks, and each time the oldest task is
   collected it starts the next module. Output order is the start order, as in the original.
2. Flag parsing. The original partitions the arguments on a leading "-", so `-j 3` would
   leave "3" as a target. This file parses the flags in order and accepts `-v`, `-j N`,
   `-jN`, `--jobs N` and `--jobs=N`.
3. The olean walk of the search path runs once, not once per target.
4. `--fresh` and the Lake manifest lookup are gone. With no target the script uses
   `StackedSVD`, which drops the `import Lake.Load.Manifest` dependency.
5. One extra final line, `KERNELREPLAY_END <n_replayed> <jobs>`, for the caller to check.

Run it from the lake project, with LEAN_PATH already set:
  lean --run scripts/KernelReplay.lean -v -j 3 [module-name-prefix ...]
-/
import Lean.CoreM
import Lean.Replay

open Lean

/-- The default size of the worker pool. Each running replay holds about 5 GB. -/
def defaultJobs : Nat := 3

/-- Unchanged from `LeanChecker.lean`: read the olean of `module`, import the oleans of its
imports, and add every constant of `module` to that environment through the kernel. -/
unsafe def replayFromImports (module : Name) : IO Unit := do
  let mFile ← findOLean module
  unless (← mFile.pathExists) do
    throw <| IO.userError s!"object file '{mFile}' of module {module} does not exist"
  -- Load all module data parts (exported, server, private)
  let mut fnames := #[mFile]
  let sFile := OLeanLevel.server.adjustFileName mFile
  if (← sFile.pathExists) then
    fnames := fnames.push sFile
    let pFile := OLeanLevel.private.adjustFileName mFile
    if (← pFile.pathExists) then
      fnames := fnames.push pFile
  let parts ← readModuleDataParts fnames
  if h : parts.size = 0 then throw <| IO.userError "failed to read module data" else
  let (mod, _) := parts[0]
  let (_, s) ← importModulesCore mod.imports |>.run
  let env ← finalizeImport s mod.imports {} 0 false false (isModule := true)
  let mut newConstants := {}
  -- Collect constants from last ("most private") part, which subsumes all prior ones
  for name in parts[parts.size-1].1.constNames, ci in parts[parts.size-1].1.constants do
    newConstants := newConstants.insert name ci
  let env' ← env.replay newConstants
  env'.freeRegions

/-- The command line: `-v`, `-j N` and the targets. Returns (verbose, jobs, targets). -/
def parseArgs (args : List String) : IO (Bool × Nat × Array String) := do
  let mut verbose := false
  let mut jobs := defaultJobs
  let mut targets : Array String := #[]
  let mut rest := args
  let readNat (s : String) : IO Nat :=
    match s.toNat? with
    | some k => pure k
    | none => throw <| IO.userError s!"KernelReplay: '{s}' is not a number of jobs"
  while !rest.isEmpty do
    match rest with
    | [] => rest := []
    | a :: as =>
      if a == "-v" || a == "--verbose" then
        verbose := true
        rest := as
      else if a == "-j" || a == "--jobs" then
        match as with
        | [] => throw <| IO.userError "KernelReplay: -j needs a number"
        | b :: bs =>
          jobs ← readNat b
          rest := bs
      else if a.startsWith "--jobs=" then
        jobs ← readNat (a.drop 7).toString
        rest := as
      else if a.startsWith "-j=" then
        jobs ← readNat (a.drop 3).toString
        rest := as
      else if a.startsWith "-j" then
        jobs ← readNat (a.drop 2).toString
        rest := as
      else if a.startsWith "-" then
        throw <| IO.userError s!"KernelReplay: unknown option '{a}'"
      else
        targets := targets.push a
        rest := as
  return (verbose, jobs, targets)

/-- Every module on the search path that one of `targets` matches. A target matches `m` when
it is a prefix of `m` or equal to `m`. The result is deduplicated and keeps discovery order. -/
def resolveTargets (targets : Array Name) : IO (Array Name) := do
  let sp ← searchPathRef.get
  let oleans ← SearchPath.findAllWithExt sp "olean"
  let mut mods : Array Name := #[]
  let mut found : Array Bool := Array.replicate targets.size false
  for path in oleans do
    if let some m := (← searchModuleNameOfFileName path sp) then
      let mut hit := false
      for i in [0:targets.size] do
        let t := targets[i]!
        if t.isPrefixOf m || t == m then
          hit := true
          found := found.set! i true
      if hit && !mods.contains m then
        mods := mods.push m
  for i in [0:targets.size] do
    unless found[i]! do
      throw <| IO.userError s!"Could not find any oleans for: {targets[i]!}"
  return mods

/-- Replay `mods`, at most `jobs` at once. Returns (no problem, number replayed). -/
unsafe def replayPool (jobs : Nat) (verbose : Bool) (mods : Array Name) : IO (Bool × Nat) := do
  let jobs := if jobs == 0 then 1 else jobs
  let n := mods.size
  let out ← IO.getStdout
  let mut tasks : Array (Task (Except IO.Error Unit)) := #[]
  let mut started := 0
  -- Fill the window.
  while started < n && started < jobs do
    tasks := tasks.push (← IO.asTask (replayFromImports mods[started]!))
    started := started + 1
  let mut done := 0
  let mut ok := true
  -- Collect the oldest task, then start one more. At most `jobs` tasks are alive.
  while done < started do
    let m := mods[done]!
    if verbose then
      out.putStrLn s!"replaying {m}"
      out.flush
    match tasks[done]!.get with
    | .error e =>
      IO.eprintln s!"leanchecker found a problem in {m}"
      IO.eprintln (toString e)
      ok := false
    | .ok _ => pure ()
    done := done + 1
    if ok && started < n then
      tasks := tasks.push (← IO.asTask (replayFromImports mods[started]!))
      started := started + 1
  return (ok, done)

/--
Replay every module that the targets match through the kernel, at most `-j N` at once.

  lean --run scripts/KernelReplay.lean -v -j 3 StackedSVD

This is not an external verifier. It detects "environment hacking", the same as `leanchecker`.
-/
unsafe def main (args : List String) : IO UInt32 := do
  initSearchPath (← findSysroot)
  let (verbose, jobs, targetStrs) ← parseArgs args
  let targetStrs := if targetStrs.isEmpty then #["StackedSVD"] else targetStrs
  let mut targets : Array Name := #[]
  for s in targetStrs do
    let mod := s.toName
    if mod.isAnonymous then
      throw <| IO.userError s!"Could not resolve module: {s}"
    targets := targets.push mod
  let mods ← resolveTargets targets
  let (ok, replayed) ← replayPool jobs verbose mods
  let out ← IO.getStdout
  out.putStrLn s!"KERNELREPLAY_END {replayed} {jobs}"
  out.flush
  return (if ok then 0 else 1)
