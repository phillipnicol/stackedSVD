#!/usr/bin/env python3
"""check_paper_edits.py - the paper edits of notes/paper_edits.md as exact, checkable replacements.

The paper (main_paper.tex, a read-only snapshot) is never edited. notes/paper_edits.md lists,
under headings of the form

    ### Exact replacement E8.1: `thm:stacksvd_weighted` (`main_paper.tex:463`)

a `Status:` line and two fenced ```latex blocks: the old text and the new text. The status is
one of:

    Status: open.
    Status: applied 2026-09-09, verbatim.
    Status: applied 2026-09-09, with changes.

An `open` replacement is still a proposal: its old text occurs verbatim (exactly once) in
main_paper.tex, and its new text is what would replace it. A replacement `applied ... verbatim`
means the user has already put the new text into main_paper.tex unchanged, so the old text no
longer occurs there. A replacement `applied ... with changes` means the user made the same
substantive change with different wording; a third fenced block, introduced by the line
"As applied (verbatim from the snapshot of DATE):", quotes what the snapshot actually says, and
a sentence after it names the difference from the proposal.

docs/THEOREMS.md quotes the amended statements under bold lead-ins of the form

    **Paper, as amended (E8).**
    **Paper, as amended (E7.1, applied 2026-09-09).**

followed by one fenced ```latex block. The label in parentheses is one or more `E<n>` group
tokens (checked against every replacement of that group) or a single `E<n>.<m>` item token
(checked against that one replacement only).

This script checks, and exits 1 on the first kind of failure it finds:
  1. every replacement has a well-formed Status line.
  2. `open`: the old block occurs exactly once in main_paper.tex, and no two open old blocks
     (or an open old block and an applied replacement's anchor, see 6) overlap; the new block
     has balanced braces and matching \\begin{...}/\\end{...} pairs.
  3. `applied ... verbatim`: the old block occurs 0 times in main_paper.tex; the new block
     occurs exactly once (this is its anchor for the overlap check).
  4. `applied ... with changes`: the old block occurs 0 times; the third ("As applied") block
     occurs exactly once (this is its anchor for the overlap check); the new block is still
     balanced, as a sanity check on the proposal text.
  5. the amended text is main_paper.tex with every OPEN replacement's new text substituted for
     its old text (applied replacements are already reflected in main_paper.tex as it stands).
     Every "Paper, as amended" block of docs/THEOREMS.md occurs verbatim in the amended text,
     and it differs from main_paper.tex itself unless every replacement named by its label is
     applied (verbatim or with changes), in which case it must occur in main_paper.tex itself.
  6. overlap: among the open old blocks and the applied replacements' anchors (an applied
     verbatim replacement's anchor is its new block's location; an applied-with-changes
     replacement's anchor is its "As applied" block's location), no open old block may overlap
     any other open old block or any applied anchor. Two applied anchors overlapping each other
     is not checked: both are read-only quotes of the same snapshot, not a pending edit, so nothing
     is at risk of being clobbered.

Usage, from the repository root:
    python3 scripts/check_paper_edits.py                 run the checks
    python3 scripts/check_paper_edits.py --apply FILE    also write the amended paper to FILE
                                                         (a scratch path; never inside the repo)
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # the repository, run from its root
PAPER = "notes/paper/main_paper.tex"
EDITS = "notes/paper_edits.md"
THEOREMS = "docs/THEOREMS.md"

HEAD = re.compile(r"^### Exact replacement (E\d+\.\d+)\b(.*)$", re.M)
FENCE = re.compile(r"```latex\n(.*?)\n```", re.S)
AMENDED = re.compile(r"^\*\*Paper, as amended \(([^)]*)\)\.\*\*[^\n]*\n\n```latex\n(.*?)\n```", re.M | re.S)
STATUS = re.compile(r"^Status: (open|applied \d{4}-\d{2}-\d{2}, verbatim|applied \d{4}-\d{2}-\d{2}, with changes)\.$", re.M)
LABEL_TOKEN = re.compile(r"E\d+(?:\.\d+)?")


def read(p):
    return open(p, encoding="utf-8").read()


def fail(msg):
    print("check_paper_edits: FAIL - " + msg)
    sys.exit(1)


def balanced(text, label):
    depth = 0
    for i, ch in enumerate(text):
        if ch == "\\" and i + 1 < len(text) and text[i + 1] in "{}":
            continue
        if i > 0 and text[i - 1] == "\\":
            continue
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth < 0:
                fail(f"{label}: a closing brace has no opening brace")
    if depth != 0:
        fail(f"{label}: braces do not balance ({depth} open)")
    stack = []
    for m in re.finditer(r"\\(begin|end)\{([^}]*)\}", text):
        if m.group(1) == "begin":
            stack.append(m.group(2))
        else:
            if not stack or stack[-1] != m.group(2):
                fail(f"{label}: \\end{{{m.group(2)}}} does not close the open environment")
            stack.pop()
    if stack:
        fail(f"{label}: environment {stack[-1]} is not closed")


def parse_status(section, tag):
    m = STATUS.search(section)
    if not m:
        fail(f"{tag}: no well-formed 'Status:' line before the first fenced block")
    label = m.group(1)
    if label == "open":
        return "open"
    if "with changes" in label:
        return "changed"
    return "verbatim"


def main():
    apply_to = None
    if len(sys.argv) == 3 and sys.argv[1] == "--apply":
        apply_to = sys.argv[2]
    elif len(sys.argv) != 1:
        print(__doc__)
        sys.exit(2)
    if not os.path.exists(PAPER) and not os.path.exists(EDITS):
        # the public release ships neither the paper snapshot nor notes/paper_edits.md
        print(f"check_paper_edits: SKIP - {PAPER} and {EDITS} are not part of this copy of the tree; "
              "the gate ran on the development tree (docs/TECHNICAL.md records the result).")
        sys.exit(0)
    paper = read(PAPER)
    edits = read(EDITS)
    heads = list(HEAD.finditer(edits))
    if not heads:
        fail(f"no '### Exact replacement' heading in {EDITS}")

    reps = {}       # tag -> dict(status, old, new, as_applied, anchor_pos, anchor_end)
    open_spans = []       # (pos, end, tag)      old-block spans of OPEN replacements
    applied_spans = []    # (pos, end, tag)      anchor spans of APPLIED replacements
    n_open = n_verb = n_changed = 0

    for k, h in enumerate(heads):
        tag = h.group(1)
        end = heads[k + 1].start() if k + 1 < len(heads) else len(edits)
        section = edits[h.end():end]
        nxt = re.search(r"^#{1,3} ", section, re.M)
        if nxt:
            section = section[:nxt.start()]
        st = parse_status(section, tag)
        blocks = FENCE.findall(section)
        expected = 3 if st == "changed" else 2
        if len(blocks) != expected:
            fail(f"{tag} (status {st}): expected {expected} ```latex blocks, found {len(blocks)}")
        old, new = blocks[0], blocks[1]
        balanced(new, tag + " new text")

        if st == "open":
            n = paper.count(old)
            if n != 1:
                fail(f"{tag} (open): old text occurs {n} times in {PAPER} (must be 1); "
                     f"first line: {old.splitlines()[0][:70]!r}")
            pos = paper.index(old)
            open_spans.append((pos, pos + len(old), tag))
            reps[tag] = dict(status=st, old=old, new=new, as_applied=None)
        elif st == "verbatim":
            n_old = paper.count(old)
            n_new = paper.count(new)
            if n_old != 0:
                fail(f"{tag} (applied verbatim): old text still occurs {n_old} times in {PAPER} "
                     f"(must be 0) - the Status line looks stale")
            if n_new != 1:
                fail(f"{tag} (applied verbatim): new text occurs {n_new} times in {PAPER} (must be 1)")
            pos = paper.index(new)
            applied_spans.append((pos, pos + len(new), tag))
            reps[tag] = dict(status=st, old=old, new=new, as_applied=None)
        else:  # changed
            as_applied = blocks[2]
            n_old = paper.count(old)
            n_ap = paper.count(as_applied)
            if n_old != 0:
                fail(f"{tag} (applied with changes): old text still occurs {n_old} times in "
                     f"{PAPER} (must be 0) - the Status line looks stale")
            if n_ap != 1:
                fail(f"{tag} (applied with changes): the 'As applied' text occurs {n_ap} times "
                     f"in {PAPER} (must be 1); first line: {as_applied.splitlines()[0][:70]!r}")
            pos = paper.index(as_applied)
            applied_spans.append((pos, pos + len(as_applied), tag))
            reps[tag] = dict(status=st, old=old, new=new, as_applied=as_applied)

        if st == "open":
            n_open += 1
        elif st == "verbatim":
            n_verb += 1
        else:
            n_changed += 1

    # overlap check: open-old spans against each other and against applied anchors;
    # applied-applied pairs are not checked (see docstring point 6).
    check_spans = sorted(open_spans + [s for s in applied_spans], key=lambda s: s[0])
    for a, b in zip(check_spans, check_spans[1:]):
        if a[1] <= b[0]:
            continue
        if a[2] in {t for t in reps if reps[t]["status"] != "open"} and \
           b[2] in {t for t in reps if reps[t]["status"] != "open"}:
            continue  # both applied: not checked
        fail(f"{a[2]} and {b[2]} overlap in {PAPER}")

    for tag in sorted(reps, key=lambda t: [int(x) for x in t[1:].split(".")]):
        r = reps[tag]
        print(f"  {tag:<6} status={r['status']:<9}"
              f" old {len(r['old'].splitlines())} lines -> new {len(r['new'].splitlines())} lines")

    # amended text: paper with every OPEN replacement's new text substituted for its old text
    amended = paper
    for pos, endp, tag in sorted(open_spans, reverse=True):
        amended = amended[:pos] + reps[tag]["new"] + amended[endp:]

    thm = read(THEOREMS)
    amended_blocks = list(AMENDED.finditer(thm))
    for m in amended_blocks:
        raw_tag, block = m.group(1), m.group(2)
        tokens = LABEL_TOKEN.findall(raw_tag)
        if not tokens:
            fail(f"THEOREMS.md amended block ({raw_tag!r}) has no E<n> label token")
        members = []
        for tok in tokens:
            if "." in tok:
                if tok not in reps:
                    fail(f"THEOREMS.md amended block ({raw_tag}): {tok} is not a known replacement")
                members.append(tok)
            else:
                group_members = [t for t in reps if t.split(".")[0] == tok]
                if not group_members:
                    fail(f"THEOREMS.md amended block ({raw_tag}): no replacement has group {tok}")
                members.extend(group_members)
        if block not in amended:
            fail(f"THEOREMS.md amended block ({raw_tag}) is not verbatim in the amended text; "
                 f"first line: {block.splitlines()[0][:70]!r}")
        all_applied = all(reps[t]["status"] != "open" for t in members)
        in_snapshot = block in paper
        if all_applied:
            if not in_snapshot:
                fail(f"THEOREMS.md amended block ({raw_tag}): every replacement of its label is "
                     f"applied, so the block must occur in {PAPER} itself; it does not")
        else:
            if in_snapshot:
                fail(f"THEOREMS.md amended block ({raw_tag}) still occurs in the unamended "
                     f"{PAPER}, but not every replacement of its label is applied")

    if apply_to:
        if os.path.abspath(apply_to).startswith(os.path.abspath(ROOT) + os.sep):
            fail("--apply target must be a scratch path outside the repository")
        open(apply_to, "w", encoding="utf-8").write(amended)
        print(f"  amended paper written to {apply_to} ({amended.count(chr(10))} lines)")

    print(f"check_paper_edits: PASS ({len(reps)} replacements: {n_verb} applied verbatim, "
          f"{n_changed} applied with changes, {n_open} open; {len(amended_blocks)} amended blocks)")


if __name__ == "__main__":
    main()
