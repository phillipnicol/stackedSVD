/-
# Layer audit driver

Checks, at the level of individual declarations rather than files, that the three Layer 2
Gaussian discharges do not use any Layer 1 theorem that consumes the hypothesis structure
they exist to prove. `RMT/Het/MPhet.lean`, `RMT/Het/Split.lean` and
`RMT/Het/Simplicity.lean` import `StackSVDWeighted.lean` at the file level, because they
need the *type* `HeteroLaw` to state what they prove. This driver confirms that the three
endpoints below never call the *proof term* of a theorem that consumes `HeteroLaw`,
`SingleTableLaw` or `HeteroEdge` (that would make Layer 2 depend on Layer 1's result, not
only on its statement). `scripts/check_layering.sh` runs this file and applies the pass or
fail rule.

Method: for each endpoint, walk `ci.type.getUsedConstants ++ ci.value?.getUsedConstants`
(the same `usedConsts` as `scripts/AxiomAudit.lean`) outward from the endpoint, expanding
only through constants declared in a `StackedSVD` module. A Mathlib or a StatsMLlib
constant cannot mention a project constant, so the walk never needs to enter one. The
result is the transitive closure of project constants that the endpoint's statement and
proof use, the endpoint included.

Run it from the lake project directory:

```
cd lean && lake env lean -j 2 ../scripts/LayerAudit.lean
```

(On the development server, `source "$LEAN_TOOLS/env.sh"` first.)

The file reads the imported environment only, so it reports the state of the oleans in
`.lake/build`, not the state of the sources (same caveat as `AxiomAudit.lean`).

## Why six imports, not the root

It imports six modules directly rather than the root `StackedSVD` file: the six reach every
`StackedSVD.*` module that the endpoints and the forbidden names need (a static import
scan of 2026-09-03: 66 project modules, including `StackSVD.lean`, `StackSVDWeighted.lean`,
`SVDStack/Main.lean`, `SVDStack/Weighted.lean`, `SVDStack/Simple.lean` and
`ThetaEst.lean`), and the root would add the rank-`r` tree and the MLE files, which no
endpoint can reach. The resolve check below fails loudly if an import goes missing.

## Two requested names do not exist verbatim

* `thm_simple_thm1_stacksvd` does not exist. `thm_simple_thm1_stacksvd_gaussian`
  (`StackSVD/Main.lean:157`) is the only theorem of that family: the paper result
  `thm:simple_thm1` has no undischarged, hypothesis-taking form for the stacksvd half in
  this tree, only the Gaussian-discharged one. It stands in below.
* `thm_simple_thm1_svdstack` does not exist either. Two Gaussian-discharged forms exist
  instead: `thm_simple_thm1_svdstack_gaussian` (`SVDStack/Main.lean:428`, stated for
  `2 ≤ M`) and `thm_simple_thm1_svdstack_gaussian_full` (`SVDStack/Simple.lean:91`, the
  paper's if/else statement for every `M ≥ 1`). Both stand in below, so a use of either
  form is caught.

`prop_dominance` (without the `_gaussian` suffix) also does not exist as a separate,
undischarged theorem; `prop:dominance` is proved only in the Gaussian-discharged form
`prop_dominance_gaussian`, which was already on the requested list.

## Output grammar

One record per line, `|` separates the fields:

```
SELFTEST <n> consumers recognized, <m> endpoints clear
LAYERAUDIT_BEGIN
ENDPOINT <full name> | closure size <k> project constants
USES <full name> | <module> | <constant>
FORBIDDEN <full name> | <forbidden constant>
CONSUMER <full name> | <constant that takes SingleTableLaw or HeteroLaw as a hypothesis>
LAYERAUDIT_END <n endpoints> <n forbidden hits>
```

A `CONSUMER` line is the list-free form of `FORBIDDEN`: a constant in the endpoint's
closure whose type has a binder of type `SingleTableLaw ...` or `HeteroLaw ...` (also under
a `∀`, as in `∀ i, (m.tbl i).SingleTableLaw (c i)`), that is a declaration that consumes one
of the two Layer 1 structures. The structures' own projections (`HeteroLaw.align` and the
like, every name with the structure name as a prefix) are exempt, since a discharge may
read a field of a structure it has built. `CONSUMER` hits count as forbidden hits in the
final line. The 13 names of `forbidden` are a fixed list; this criterion catches a consumer
the list does not name. The `SELFTEST` line comes first: the criterion is run on five known
consumers (must be `true`) and on the three endpoints (must be `false`) before the audit,
and a failure is a Lean error, so a broken criterion cannot pass the gate silently.

A `USES` line reports a constant in the endpoint's closure whose defining module is
`StackedSVD.StackSVDWeighted`, `StackedSVD.StackSVD`, or a module under
`StackedSVD.SVDStack`: the Layer 1 file *definitions* (structures, `def`s, such as
`HeteroLaw`, `stackGramW`, `stackPerfW`) that a Gaussian discharge is allowed to see, as
opposed to the *theorems* named in `forbidden`, which it must not call. `USES` lines for
one endpoint are sorted by constant name and free of duplicates (the closure is a set, so
a constant appears at most once).

A `FORBIDDEN` line reports a forbidden name found in an endpoint's closure; none are
expected. The driver throws an error, instead of silently printing a wrong count, when an
endpoint or a forbidden name fails to resolve in the imported environment.
-/
import StackedSVD.RMT.Full
import StackedSVD.RMT.Het.Sup
import StackedSVD.RMT.Het.EdgeSharp
import StackedSVD.SVDStack.Weighted
import StackedSVD.ThetaEst
import StackedSVD.SVDStack.Simple

open Lean Elab Command

namespace LayerAudit

/-- Module prefix of the project. -/
def projPrefix : Name := `StackedSVD

/-- Every constant that occurs in the type or in the value (the proof term). Identical to
`AxiomAudit.usedConsts`. -/
def usedConsts (ci : ConstantInfo) : Array Name :=
  ci.type.getUsedConstants ++
    (match ci.value? (allowOpaque := true) with
     | some v => v.getUsedConstants
     | none => #[])

/-- Transitive closure of the project constants that `n` uses, `n` included. The walk
stops at every constant outside `proj`, since a Mathlib or a StatsMLlib constant cannot
mention a project constant. Identical in shape to `AxiomAudit.projDeps`. -/
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

/-- Module that declares `n`. Built from `Environment.getModuleIdxFor?` the same way
`ImportGraph.Lean.Environment.getModuleFor?` does; that helper is not imported here, to
keep this driver self-contained (one extra project import for one lookup). `none` only for
a constant outside every imported module, which does not happen for a name already found
by `env.find?`. -/
def moduleOf (env : Environment) (n : Name) : Option Name :=
  (env.getModuleIdxFor? n).bind fun idx => env.header.moduleNames[idx.toNat]?

/-- `true` for the module of a Layer 1 file: `StackSVDWeighted.lean`, `StackSVD.lean`, or
anything under `SVDStack/`. These are the files whose *definitions* a Layer 2 discharge is
allowed to use. -/
def isLayer1Module (m : Name) : Bool :=
  m == `StackedSVD.StackSVDWeighted || m == `StackedSVD.StackSVD ||
    m == `StackedSVD.StackSVD.Main || m == `StackedSVD.StackSVD.Weighted ||
    m == `StackedSVD.RankR.Unweighted || m == `StackedSVD.RankR.WeightedMain ||
    m == `StackedSVD.RankR.SubspaceMain ||
    (`StackedSVD.SVDStack).isPrefixOf m

/-- The Layer 2 targets under audit: the Gaussian discharge of `SingleTableLaw`,
`HeteroLaw` and `HeteroEdge`. -/
def endpoints : List Name :=
  [ `StackedSVD.SpikedModel.singleTableLaw_of_gaussian,
    `StackedSVD.MultiTableModel.heteroLaw_of_gaussian,
    `StackedSVD.MultiTableModel.heteroEdge_of_gaussian ]

/-- Layer 1 theorems that consume one of the three hypothesis structures to prove the
paper's probabilistic claims. None may appear in an endpoint's closure. See the file doc
comment above for the two names substituted for a non-existent request. -/
def forbidden : List Name :=
  [ `StackedSVD.MultiTableModel.prop_stacksvd_general,
    `StackedSVD.MultiTableModel.prop_stacksvd_general_gaussian,
    `StackedSVD.MultiTableModel.thm_stacksvd_weighted,
    `StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian,
    `StackedSVD.MultiTableModel.thm_svd_stack_general,
    `StackedSVD.MultiTableModel.thm_svd_stack_general_gaussian,
    `StackedSVD.MultiTableModel.thm_svdstack_weighted,
    `StackedSVD.MultiTableModel.thm_svdstack_weighted_gaussian,
    `StackedSVD.MultiTableModel.prop_dominance_gaussian,
    `StackedSVD.MultiTableModel.thm_theta_est,
    `StackedSVD.MultiTableModel.thm_simple_thm1_stacksvd_gaussian,
    `StackedSVD.MultiTableModel.thm_simple_thm1_svdstack_gaussian,
    `StackedSVD.MultiTableModel.thm_simple_thm1_svdstack_gaussian_full ]

/-- The two Layer 1 hypothesis structures. `HeteroEdge` (`RMT/Het/R4het.lean`) is an
interface inside Layer 2 and is not one of them. -/
def structs : List Name :=
  [ `StackedSVD.SpikedModel.SingleTableLaw,
    `StackedSVD.MultiTableModel.HeteroLaw ]

/-- Head constant of an expression after its leading `∀` binders are stripped, so that a
hypothesis `∀ i, (m.tbl i).SingleTableLaw (c i)` is seen as `SingleTableLaw`. Loose bound
variables are harmless here: only the head is read. -/
partial def headAfterForalls : Expr → Option Name
  | .forallE _ _ b _ => headAfterForalls b
  | e => match e.getAppFn with
    | .const n _ => some n
    | _ => none

/-- `true` when some binder of the type `ty` has one of `structs` as its type, that is when
the declaration takes a Layer 1 structure as a hypothesis. A declaration that only
*concludes* the structure (the endpoints themselves) is not a consumer. -/
partial def consumesStruct : Expr → Bool
  | .forallE _ t b _ =>
      (match headAfterForalls t with
       | some n => structs.contains n
       | none => false) || consumesStruct b
  | _ => false

/-- The structure's own projections, recursors and constructor are exempt from the
consumer criterion: their names have the structure name as a prefix. -/
def isStructInternal (n : Name) : Bool :=
  structs.any fun s => s.isPrefixOf n

def run : CommandElabM Unit := do
  let env ← getEnv
  let header := env.header
  let mods := (header.moduleNames.zip header.moduleData).filter fun (mn, _) =>
    projPrefix.isPrefixOf mn
  let projConsts : NameSet := mods.foldl (init := {}) fun s (_, md) =>
    md.constNames.foldl (init := s) fun s c => s.insert c
  -- Resolve check: every endpoint, every forbidden name and both structures must exist
  -- in this environment, or the audit below would silently under-report.
  let mut unresolved : Array Name := #[]
  for n in endpoints ++ forbidden ++ structs do
    if (env.find? n).isNone then
      unresolved := unresolved.push n
  unless unresolved.isEmpty do
    throwError "LayerAudit: unresolved name(s), fix the import list or the spelling: {unresolved}"
  -- Positive control for the consumer criterion: five Layer 1 theorems that take a
  -- `SingleTableLaw` or a `HeteroLaw` hypothesis must register as consumers, and the
  -- three endpoints, which conclude a structure but take none, must not. A failure
  -- here is a bug in `consumesStruct`, so it is a Lean error (exit 2 in the gate).
  let positives : List Name :=
    [ `StackedSVD.MultiTableModel.prop_stacksvd_general,
      `StackedSVD.MultiTableModel.thm_stacksvd_weighted,
      `StackedSVD.MultiTableModel.thm_svd_stack_general,
      `StackedSVD.MultiTableModel.thm_svdstack_weighted,
      `StackedSVD.MultiTableModel.thm_theta_est ]
  for n in positives do
    match env.find? n with
    | some ci =>
        unless consumesStruct ci.type do
          throwError "LayerAudit: self-test failed, {n} is not seen as a consumer"
    | none => throwError "LayerAudit: self-test name {n} does not resolve"
  for n in endpoints do
    match env.find? n with
    | some ci =>
        if consumesStruct ci.type then
          throwError "LayerAudit: self-test failed, endpoint {n} is seen as a consumer"
    | none => pure ()
  IO.println s!"SELFTEST {positives.length} consumers recognized, {endpoints.length} endpoints clear"
  IO.println "LAYERAUDIT_BEGIN"
  let mut totalForbidden := 0
  for ep in endpoints do
    let closure := projDeps env projConsts ep {}
    IO.println s!"ENDPOINT {ep} | closure size {closure.size} project constants"
    let mut uses : Array (Name × Name) := #[]
    for c in closure do
      match moduleOf env c with
      | some m => if isLayer1Module m then uses := uses.push (m, c)
      | none => pure ()
    let sorted := uses.qsort fun a b => (compare (toString a.2) (toString b.2)).isLT
    for (m, c) in sorted do
      IO.println s!"USES {ep} | {m} | {c}"
    for f in forbidden do
      if closure.contains f then
        totalForbidden := totalForbidden + 1
        IO.println s!"FORBIDDEN {ep} | {f}"
    let mut consumers : Array Name := #[]
    for c in closure do
      if c != ep && !(isStructInternal c) then
        match env.find? c with
        | some ci => if consumesStruct ci.type then consumers := consumers.push c
        | none => pure ()
    for c in consumers.qsort (fun a b => (compare (toString a) (toString b)).isLT) do
      totalForbidden := totalForbidden + 1
      IO.println s!"CONSUMER {ep} | {c}"
  IO.println s!"LAYERAUDIT_END {endpoints.length} {totalForbidden}"

end LayerAudit

run_cmd do LayerAudit.run
