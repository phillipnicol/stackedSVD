/-
Copyright (c) 2026 Tavor Baharav. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tavor Baharav, Claude
-/
import StackedSVD.SVDStack.Defs
import StackedSVD.SVDStack.Deterministic
import StackedSVD.SVDStack.Gram
import StackedSVD.SVDStack.Main
import StackedSVD.SVDStack.Weighted
import StackedSVD.SVDStack.Simple
import StackedSVD.SVDStack.Rayleigh
import StackedSVD.SVDStack.EntrywiseSigned
import StackedSVD.SVDStack.Inad

/-!
# `thm:svd_stack_general`: unweighted svdstack

This file is the entry point of the directory `StackedSVD/SVDStack/`. It imports the
modules of the directory and adds nothing, so `import StackedSVD.SVDStack` keeps its old
meaning. Added 2026-09-02: `EntrywiseSigned` (the paper's signed coordinate form of
`lem:entrywise_conv_eigenvec`, item L2 of the external statement audit) and `Inad` (the
svdstack clauses of `prop:binarystacksvd_inadmissable` on the model, item L3).

STATUS 2026-08-29 (evening): every result of the directory is proved; no deferred proof is
left. `topSimple_svdstackGram` was restated with probability tending to one (user decision 5
of `notes/archive/thm_svd_stack_general.md`; the a.s. form needs a density for the joint law of the
top eigenvectors, which no hypothesis gives). `thm_svd_stack_general` is assembled from the
general lemmas of the section `EntrywiseEigenvec` (`topProj_overlap_tendsto`,
`lamMax_tendstoInProb`, `topSimple_whp_of_tendsto`, `lamMax_gt_half_whp_of_tendsto`), which
take any symmetric limit matrix so that the weighted theorem can reuse them, and the paper's
form `|⟨v̂_svdstack, v⟩|² → ...` is `thm_svd_stack_general_inner`
(`notes/archive/agent_reports/proof_svdstack_final.md`).

svdstack takes the top right singular vector `v̂_i` of each table and returns the top right
singular vector of the `M × d` matrix `Ṽ` whose rows are the `v̂_iᵀ`. Its limit is
`(βᵀ v_max(A_β))² / λ_max(A_β)` with `A_β = β βᵀ + diag(1 - β_i²)`.

## The four modules

1. `StackedSVD.SVDStack.Defs`: `beta`, `Abeta`, `vMax`, the spectral helpers, the closed form
   `svdstackPerfClosed` with `svdstackPerf_eq_closed`, and the per-table and svdstack objects
   `tableGram`, `vhat`, `P`, `Vt`, `svdstackGram`, `svdstackEst`, `svdstackPerf`.
2. `StackedSVD.SVDStack.Deterministic`: `beta_mem_Ico`, `one_le_lamMax_Abeta`, `abeta_gap`,
   and the limit lemmas of the section `EntrywiseEigenvec` for any symmetric limit matrix.
3. `StackedSVD.SVDStack.Gram`: `svdstackGram_eq`, the `delocDir` machinery,
   `lem_delocalization`, `gram`, `align`, `goodEvent` and `topSimple_svdstackGram`.
4. `StackedSVD.SVDStack.Main`: `thm_svd_stack_general`, `thm_svd_stack_general_inner` and
   `thm_svd_stack_general_zero`.

All limit statements about `Ṽ` are signed, not squared. A squared `gram` does not determine
the conclusion: two limit matrices with the same squared entries give svdstack limits
`0.529412` and `0` (`notes/archive/audit_scope_2026-08-29.md` section 6.1). The sign convention of
`vhat` is what makes the signed form correct.
-/
