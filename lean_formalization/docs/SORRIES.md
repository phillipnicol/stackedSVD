# docs/SORRIES.md

Every `sorry` under `lean/StackedSVD/` has one row here. Keep in sync with `grep -rn sorry lean/StackedSVD/`.

`scripts/check_sorries.sh` is the gate. It compares this table with the tree after it blanks
every comment, so a `sorry` named inside a docstring is not a site.

The table is empty. The last open proof was the one combinatorial input of Stage 3 of
`notes/stage3_edge.md` (the sharp upper edge of the noise Gram matrix at a general law,
route C, the moment method with a truncation): the Furedi-Komlos count,
`cellCard_mul_le_pos` (a closed walk of length `2k` on the complete bipartite graph with no
entry visited exactly once, and with one vertex fewer than a tree walk, costs at most
`(2k)^12 / min(n, d)` per missing vertex). It is now proved, in three modules of
`RMT/General/Edge/`: `Code.lean` (the six-part step code of a walk and its injectivity),
`Dyck.lean` (the Dyck-word tree walks and the Narayana lower bound), and `CountBound.lean`
(the three arithmetic lemmas and the assembly `cellCard_mul_le_pos`).

The count also has an independent numeric witness. Exact enumeration checks it for `k` up to
7 and for `(n, d)` from `(k+1, k+1)` to `(10^12, 10^9)`: 0 violations, and it still holds
with `(2k)^2` in place of `(2k)^12`, so the exponent has a margin of `(2k)^10`. Section 6 of
`notes/STAGE3_CAMPAIGN.md` names a fallback route, the typed-edge count of Yin, Bai and
Krishnaiah 1988; it was not needed.

| File | Declaration | Reason | Since |
|------|-------------|--------|-------|
