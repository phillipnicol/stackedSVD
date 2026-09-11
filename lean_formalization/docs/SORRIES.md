# docs/SORRIES.md

Every `sorry` under `lean/StackedSVD/` has one row here. Keep in sync with `grep -rn sorry lean/StackedSVD/`.

`scripts/check_sorries.sh` is the gate. It compares this table with the tree after it blanks
every comment, so a `sorry` named inside a docstring is not a site.

The table is empty since 2026-09-10 evening. Its last row was the one combinatorial input of
Stage 3 of `notes/stage3_edge.md` (the sharp upper edge of the noise Gram matrix at a general
law, route C, the moment method with a truncation): the Furedi-Komlos count,
`cellCard_mul_le_pos` (a closed walk of length `2k` on the complete bipartite graph with no
entry visited exactly once, and with one vertex fewer than a tree walk, costs at most
`(2k)^12 / min(n, d)` per missing vertex). It was stated and used, so every other Stage 3
result was proved modulo it until it too was discharged. The count is checked by exact
enumeration for `k` up to 7 and for `(n, d)` from `(k+1, k+1)` to `(10^12, 10^9)`: 0
violations, and it still holds with `(2k)^2` in place of `(2k)^12`, so the exponent has a
margin of `(2k)^10`. Risk 1 of `notes/STAGE3_CAMPAIGN.md` (the Stage 3 campaign note)
section 6 names the fallback route, the typed-edge count of Yin, Bai and Krishnaiah 1988; it
was not needed. The row was the last open item of the 34 Stage 3 statements, added
2026-09-10 by workflow step 3 of `CLAUDE.md` on the branch `stage3` (the statements OK by
trust, the user's message of 2026-09-10, with a Fable audit of every statement recorded in
the campaign note) and discharged the same day: 33 of 34, each by the unit of the campaign
note named in its row: unit C, the comparison of the truncated trace with the Gaussian one
(`Edge/Compare.lean`), and unit T, the truncation (`Edge/Trunc.lean`), in the late morning;
unit L, the discarded part (`Edge/Sparse.lean`), unit E, the trace expansion over closed
walks (`Edge/Trace.lean`), unit G, the Gaussian trace bound (`Edge/Gaussian.lean`), and unit
K, the arithmetic of the moment order and the Markov step (`Edge/Arith.lean`,
`Edge/Markov.lean`), after them; and unit X, the excess count (`Edge/Excess.lean`). The 34th,
the Furedi-Komlos count itself, was proved last, 2026-09-10 evening, by the X8 plan
(`notes/x8_plan.md`, the proof plan for the count) in three new modules of
`RMT/General/Edge/`: `Code.lean` (the six-part step code of a walk and its injectivity),
`Dyck.lean` (the Dyck-word tree walks and the Narayana lower bound), and `CountBound.lean`
(the three arithmetic lemmas and the assembly `cellCard_mul_le_pos`); the exponent of the
count changed from 6 to 12 in the process. Before them the table was empty (2026-09-09
evening). Its
last rows before that were the 12 Stage 1 statements of
`notes/archive/prop_single_table_general.md` (the non-Gaussian single-table law, user OK of the
statements 2026-09-09), added by workflow step 3 of `CLAUDE.md` and discharged the same day,
each by the unit of section 5 of that note named in its row: the item S row
(`singleTableLaw_topSimple_of_general`, unit G8) in the morning; the two endpoint rows
(`resolventLimits_of_general`, `singleTableLaw_of_general`, unit G10, moved to
`RMT/General/Sup.lean`) and the six rows of `General/Layer1.lean` (unit G11) in the evening,
each proved against the three stubs that then remained; the six complex forms
(`resolventFormsC_of_general`, unit F3, now `resolventFormsC_of_general'` in
`RMT/General/FormsGeneral.lean`) and the subcritical align field
(`align_tendstoInProb_of_subcritical_general`, unit D1c, now
`align_tendstoInProb_of_subcritical_general'` in `RMT/General/DelocAlign.lean`, with the
three added hypotheses `RL`, `[SigmaFinite ν]` and `hac`) later that evening; and last the
delocalization field (`delocUniform_of_general`, unit G7, now proved under the same name in
`RMT/General/Sup.lean` from `delocUniform_of_lamMax` in `RMT/General/DelocUniform.lean`).
Before them the table was empty; the rows before those were the 26 statements of item F18b
(`notes/archive/F18b_plan.md`, the Gaussian discharge of the `hlaw` hypothesis of
`prop_singleweight_suboptimality_of_law`), added 2026-09-06 by workflow step 3 of `CLAUDE.md`
after the user's OK of the statements on 2026-09-06 and discharged on 2026-09-06 and
2026-09-07. Before them the table was empty.
| File | Declaration | Reason | Since |
|------|-------------|--------|-------|
