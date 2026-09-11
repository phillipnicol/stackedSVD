/-
# Axiom audit driver

Prints the axiom dependencies of every declaration of the `StackedSVD` modules,
private declarations included. `scripts/check_axioms.sh` parses the output and
applies the pass or fail rule.

Run it from the lake project directory:

```
cd lean && lake env lean -j 3 ../scripts/AxiomAudit.lean
```

(On the development server, `source "$LEAN_TOOLS/env.sh"` first.)

The file reads the imported environment only, so it reports the state of the oleans
in `.lake/build`, not the state of the sources. See `scripts/README.md`.

Output grammar (one record per line, `|` separates the fields):

```
AXIOMAUDIT_BEGIN
MODULE <module name>
DECL <declaration> | <axiom>,<axiom>,...        -- `-` when there is no axiom
SORRYDEPS <declaration> | <name>,<name>,...     -- only after a `DECL` with `sorryAx`
AXIOMAUDIT_END <reported> <skipped>
```

A private declaration is printed under its user facing name, with the suffix
` (private)` after the name. `SORRYDEPS` lists the project constants that carry a
`sorry` in their own body and that the declaration uses. It lets the shell tell a
tracked `sorry` from a new one.
-/
import StackedSVD
import Lean.Util.CollectAxioms

open Lean Elab Command

namespace AxiomAudit

/-- Module prefix of the project. -/
def projPrefix : Name := `StackedSVD

/-- Vendored third party code. Skipped unless `AXIOM_AUDIT_VENDOR` is set to `1`. -/
def vendorPrefix : Name := `StackedSVD.Vendor

/-- Last name components that the compiler generates. -/
def noiseComponents : List String :=
  ["rec", "recOn", "casesOn", "brecOn", "below", "ibelow", "binductionOn", "ndrec",
   "noConfusion", "noConfusionType", "injEq", "sizeOf_spec", "sizeOf_inst", "toCtorIdx",
   "eq_def", "sunfold", "induct", "fun_cases", "ofNat", "mk", "inj", "ctorIdx"]

/-- True for an auxiliary or compiler generated name. `Name.isInternalDetail` covers
`_auxLemma`, `.proof_1`, `.match_1`, `.eq_1` and every `_`-prefixed part. Apply it to the
user facing name, since a private name always starts with `_private`. -/
def isNoise (n : Name) : Bool :=
  if n.isInternalDetail || n.hasMacroScopes then true
  else match n with
    | .str _ s => noiseComponents.contains s
    | _ => true

/-- Head symbol of a type, after the binders. -/
partial def resultHead : Expr → Name
  | .forallE _ _ b _ => resultHead b
  | e => e.getAppFn.constName?.getD Name.anonymous

/-- `Decidable` instances are generated noise for this audit. -/
def isDecidableLike (ci : ConstantInfo) : Bool :=
  let h := resultHead ci.type
  h == ``Decidable || h == ``DecidableEq || h == ``DecidablePred

/-- Every constant that occurs in the type or in the value. -/
def usedConsts (ci : ConstantInfo) : Array Name :=
  ci.type.getUsedConstants ++
    (match ci.value? (allowOpaque := true) with
     | some v => v.getUsedConstants
     | none => #[])

/-- True when the body of `n` mentions `sorryAx` itself, so the `sorry` is in this
declaration and not in something that it uses. -/
def hasOwnSorry (env : Environment) (n : Name) : Bool :=
  match env.find? n with
  | some ci => (usedConsts ci).contains ``sorryAx
  | none => false

/-- Transitive closure of the project constants that `n` uses, `n` included.
The walk stops at every constant outside `proj`. -/
partial def projDeps (env : Environment) (proj : NameSet) (n : Name) (acc : NameSet) :
    NameSet :=
  if acc.contains n then acc
  else
    let acc := acc.insert n
    match env.find? n with
    | none => acc
    | some ci =>
      (usedConsts ci).foldl (init := acc) fun a m =>
        if proj.contains m then projDeps env proj m a else a

/-- The project declarations that carry the `sorry` on which `n` depends. -/
def sorrySources (env : Environment) (proj : NameSet) (n : Name) : Array Name := Id.run do
  let mut out : Array Name := #[]
  for m in projDeps env proj n {} do
    if hasOwnSorry env m then out := out.push (privateToUserName m)
  return out

def fmtNames (a : Array Name) : String :=
  if a.isEmpty then "-" else String.intercalate "," (a.toList.map toString)

def run (withVendor : Bool) : CommandElabM Unit := do
  let env ← getEnv
  let header := env.header
  let mods := (header.moduleNames.zip header.moduleData).filter fun (mn, _) =>
    projPrefix.isPrefixOf mn && (withVendor || !vendorPrefix.isPrefixOf mn)
  -- Every constant declared in a project module, private names included. `projDeps`
  -- uses this instead of a name prefix, because a private name reads `_private....`.
  let projConsts : NameSet := mods.foldl (init := {}) fun s (_, md) =>
    md.constNames.foldl (init := s) fun s c => s.insert c
  let mut reported := 0
  let mut skipped := 0
  IO.println "AXIOMAUDIT_BEGIN"
  for (mname, mdata) in mods do
    IO.println s!"MODULE {mname}"
    for c in mdata.constNames do
      let user := privateToUserName c
      let mark := if isPrivateName c then " (private)" else ""
      if isNoise user then
        skipped := skipped + 1
        continue
      match env.find? c with
      | none => skipped := skipped + 1
      | some ci =>
        if isDecidableLike ci then
          skipped := skipped + 1
        else
          let axs ← collectAxioms c
          reported := reported + 1
          IO.println s!"DECL {user}{mark} | {fmtNames axs}"
          if axs.contains ``sorryAx then
            IO.println s!"SORRYDEPS {user}{mark} | {fmtNames (sorrySources env projConsts c)}"
  IO.println s!"AXIOMAUDIT_END {reported} {skipped}"

end AxiomAudit

run_cmd do
  let withVendor := (← IO.getEnv "AXIOM_AUDIT_VENDOR") == some "1"
  AxiomAudit.run withVendor
