/-
AuditSignatures.lean - elaborated signatures and axiom dependencies for AUDIT_PACK_V3 and V4.

This file is NOT part of the build graph. It exists to answer two requests of the external
audit of 2026-08-30 (`notes/archive/external_audit_2026-08-30.md`, action items 14 and 17):

  14. show the ambient variables and instances that a source snippet inside a `section`
      hides, by printing the elaborated constant instead of the source text;
  17. show `#print axioms` for every headline theorem.

Run it from the lake project directory. The import resolves through `lake env`, so the file
does not need to sit under the package root:

  cd lean && lake env lean -j 3 ../scripts/AuditSignatures.lean

(On the development server, `source "$LEAN_TOOLS/env.sh"` first.)

Every message is self-identifying: `#check @f` prints `@f : ...` and `#print axioms f`
prints `'f' depends on axioms: [...]`.
-/
import StackedSVD

set_option linter.style.header false

/-! ## Part 1: elaborated signatures (`#check @f`) -/

/-! ### The three hypothesis structures -/

#check @StackedSVD.SpikedModel.SingleTableLaw
#check @StackedSVD.MultiTableModel.HeteroLaw
#check @StackedSVD.MultiTableModel.HeteroEdge

/-! ### Headline theorems, Gaussian noise -/

#check @StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian
#check @StackedSVD.MultiTableModel.heteroLaw_of_gaussian
#check @StackedSVD.MultiTableModel.heteroEdge_of_gaussian
#check @StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian_margin
#check @StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian_margin_inner
#check @StackedSVD.MultiTableModel.thm_svdstack_weighted_gaussian
#check @StackedSVD.MultiTableModel.thm_svd_stack_general_gaussian
#check @StackedSVD.MultiTableModel.prop_stacksvd_general_gaussian
#check @StackedSVD.SpikedModel.singleTableLaw_of_gaussian
#check @StackedSVD.MultiTableModel.thm_theta_est_gaussian
#check @StackedSVD.MultiTableModel.lem_delocalization

/-! ### Layer 1 theorems (a law structure is a hypothesis) -/

#check @StackedSVD.MultiTableModel.prop_stacksvd_general
#check @StackedSVD.MultiTableModel.thm_stacksvd_weighted
#check @StackedSVD.MultiTableModel.thm_stacksvd_weighted_general
#check @StackedSVD.MultiTableModel.thm_svd_stack_general
#check @StackedSVD.MultiTableModel.thm_svdstack_weighted
#check @StackedSVD.MultiTableModel.thm_theta_est
#check @StackedSVD.MultiTableModel.stackPerfW_binary_inad
#check @StackedSVD.MultiTableModel.stackPerfW_opt_inad

/-! ### The scalar declarations that the word "optimal" rests on -/

#check @StackedSVD.Scalars.L_le_opt
#check @StackedSVD.Scalars.L_optW_eq
#check @StackedSVD.Scalars.assumption4_optW_iff

/-! ### Campaign E: the exact heteroscedastic edge

`MPhet.edgeObjective` and its scalar bridge lemmas are pure `ℝ` statements with no ambient
model context, so their source quotes in sections 4.7 and 5.7b are already complete. The same
holds for the `Scalars` declarations of section 6; only the three that the word "optimal"
rests on are checked below. -/

#check @StackedSVD.integral_opNorm_mul_le_sharp
#check @StackedSVD.tendsto_measure_lamMax_le_of_bound

/-! ### The rank-one Gaussian facade (commit 6705361)

Every paper claim of sections 3 to 6 has a named Gaussian endpoint whose hypotheses are model
hypotheses only. The 17 declarations below are those endpoints and the two bundles.

Section 7 and Appendix E had no entry here when this file was written, because the rank-`r`
work sat in `parked/RankR/`, outside the package. That is no longer true: `parked/` is gone,
every rank-`r` module is in the build, and Parts 2 to 8 below cover it. -/

/-! #### stackSVD, unweighted (`StackSVD.lean`) -/

#check @StackedSVD.MultiTableModel.thm_simple_thm1_stacksvd_gaussian

/-! #### SVDstack, unweighted (`SVDStack/Main.lean`) -/

#check @StackedSVD.MultiTableModel.thm_svd_stack_general_zero_gaussian
#check @StackedSVD.MultiTableModel.thm_svd_stack_general_inner_gaussian
#check @StackedSVD.MultiTableModel.thm_simple_thm1_svdstack_gaussian
#check @StackedSVD.MultiTableModel.thm_simple_thm1_svdstack_gaussian_full

/-! #### SVDstack, weighted (`SVDStack/Weighted.lean`) -/

#check @StackedSVD.MultiTableModel.thm_svdstack_weighted_general_gaussian
#check @StackedSVD.MultiTableModel.thm_svdstack_weighted_paper_gaussian
#check @StackedSVD.MultiTableModel.thm_svdstack_weighted_inner_gaussian
#check @StackedSVD.MultiTableModel.thm_svdstack_weighted_zero_gaussian
#check @StackedSVD.MultiTableModel.thm_svdstack_weighted_gaussian_opt

/-! #### stackSVD, weighted (`RMT/Het/Sup.lean`, `StackSVDWeighted.lean`) -/

#check @StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian_inner
#check @StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian_smul
#check @StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian_opt
#check @StackedSVD.MultiTableModel.prop_dominance_gaussian
#check @StackedSVD.MultiTableModel.stackPerfW_binary_inad_gaussian
#check @StackedSVD.MultiTableModel.stackPerfW_opt_inad_gaussian
#check @StackedSVD.MultiTableModel.thm_stacksvd_binary_optimal_svd_stack_gaussian

/-! #### Strict comparisons (`StrictFacades.lean`, after the third external review, P2) -/

#check @StackedSVD.MultiTableModel.thm_stacksvd_binary_optimal_svd_stack_gaussian_strict
#check @StackedSVD.MultiTableModel.prop_dominance_gaussian_strict
#check @StackedSVD.MultiTableModel.prop_dominance_gaussian_strict_svdstack
#check @StackedSVD.MultiTableModel.prop_dominance_gaussian_strict_unweighted

/-! #### Uniform Rayleigh bound (`SVDStack/Rayleigh.lean`, P1): no weight beats `S/(S+1)` -/

#check @StackedSVD.MultiTableModel.rowBound
#check @StackedSVD.MultiTableModel.svdstackPerfW_le_rowBound
#check @StackedSVD.MultiTableModel.svdstackPerfW_uniform_bound
#check @StackedSVD.MultiTableModel.svdstackPerfW_le_opt_whp
#check @StackedSVD.MultiTableModel.thm_svdstack_weighted_gaussian_opt_full
#check @StackedSVD.svdstackLimitOpt_eq_of_single

/-! #### `app:wstacksvd_mle` (`MLE.lean`, D5): the marginal MLE is weighted stackSVD -/

#check @StackedSVD.MultiTableModel.mleLogLik
#check @StackedSVD.MultiTableModel.mleLogLik_eq
#check @StackedSVD.MultiTableModel.mleLogLik_le_iff
#check @StackedSVD.MultiTableModel.mleLogLik_le_of_mem_topSpace
#check @StackedSVD.inner_toOp_self_eq_lamMax_iff
#check @StackedSVD.MultiTableModel.mem_topSpace_of_mleLogLik_max
#check @StackedSVD.MultiTableModel.mleLogLik_max_iff_mem_topSpace

/-! #### `remark:stack_outperform_svd`, `remark:svd_outperform_stack` (`Remarks.lean`, P4) -/

#check @StackedSVD.svdstackLimit_pair
#check @StackedSVD.remark_three_stack_tendsto_zero
#check @StackedSVD.MultiTableModel.remark_stack_outperform_svd_stack
#check @StackedSVD.MultiTableModel.remark_stack_outperform_svd_svdstack
#check @StackedSVD.MultiTableModel.remark_svd_outperform_stack_pair_stack
#check @StackedSVD.MultiTableModel.remark_svd_outperform_stack_pair_svdstack
#check @StackedSVD.MultiTableModel.remark_svd_outperform_stack_three_svdstack
#check @StackedSVD.MultiTableModel.remark_svd_outperform_stack_three_stack
#check @StackedSVD.MultiTableModel.remark_svd_outperform_stack_three_binary
#check @StackedSVD.MultiTableModel.remark_svd_outperform_stack_two
#check @StackedSVD.svdstackLimitOpt_zero
#check @StackedSVD.MultiTableModel.remark_stack_outperform_svd_svdstack_uniform
#check @StackedSVD.MultiTableModel.remark_svd_outperform_stack_two_svdstack_unweighted

/-! ### `lem:secular_equation`: the multiplicity clause and the block bridge -/

#check @StackedSVD.Secular.card_sub_one_le_finrank_eigenspace
#check @StackedSVD.Secular.secularDiag_stack

/-! ## Part 2: axiom dependencies (`#print axioms f`) -/

#print axioms StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian
#print axioms StackedSVD.MultiTableModel.heteroLaw_of_gaussian
#print axioms StackedSVD.MultiTableModel.heteroEdge_of_gaussian
#print axioms StackedSVD.MultiTableModel.thm_svdstack_weighted_gaussian
#print axioms StackedSVD.MultiTableModel.prop_stacksvd_general_gaussian
#print axioms StackedSVD.MultiTableModel.thm_svd_stack_general_gaussian
#print axioms StackedSVD.SpikedModel.singleTableLaw_of_gaussian
#print axioms StackedSVD.MultiTableModel.thm_theta_est_gaussian
#print axioms StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian_margin
#print axioms StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian_margin_inner
#print axioms StackedSVD.MultiTableModel.lem_delocalization
#print axioms StackedSVD.Scalars.L_le_opt
#print axioms StackedSVD.Scalars.L_optW_eq
#print axioms StackedSVD.Scalars.assumption4_optW_iff
#print axioms StackedSVD.integral_opNorm_mul_le_sharp
#print axioms StackedSVD.tendsto_measure_lamMax_le_of_bound

/-! ### The rank-one Gaussian facade -/

#print axioms StackedSVD.MultiTableModel.thm_simple_thm1_stacksvd_gaussian
#print axioms StackedSVD.MultiTableModel.thm_svd_stack_general_zero_gaussian
#print axioms StackedSVD.MultiTableModel.thm_svd_stack_general_inner_gaussian
#print axioms StackedSVD.MultiTableModel.thm_simple_thm1_svdstack_gaussian
#print axioms StackedSVD.MultiTableModel.thm_simple_thm1_svdstack_gaussian_full
#print axioms StackedSVD.MultiTableModel.thm_svdstack_weighted_general_gaussian
#print axioms StackedSVD.MultiTableModel.thm_svdstack_weighted_paper_gaussian
#print axioms StackedSVD.MultiTableModel.thm_svdstack_weighted_inner_gaussian
#print axioms StackedSVD.MultiTableModel.thm_svdstack_weighted_zero_gaussian
#print axioms StackedSVD.MultiTableModel.thm_svdstack_weighted_gaussian_opt
#print axioms StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian_inner
#print axioms StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian_smul
#print axioms StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian_opt
#print axioms StackedSVD.MultiTableModel.prop_dominance_gaussian
#print axioms StackedSVD.MultiTableModel.stackPerfW_binary_inad_gaussian
#print axioms StackedSVD.MultiTableModel.stackPerfW_opt_inad_gaussian
#print axioms StackedSVD.MultiTableModel.thm_stacksvd_binary_optimal_svd_stack_gaussian
#print axioms StackedSVD.MultiTableModel.thm_stacksvd_binary_optimal_svd_stack_gaussian_strict
#print axioms StackedSVD.MultiTableModel.prop_dominance_gaussian_strict
#print axioms StackedSVD.MultiTableModel.prop_dominance_gaussian_strict_svdstack
#print axioms StackedSVD.MultiTableModel.prop_dominance_gaussian_strict_unweighted
#print axioms StackedSVD.MultiTableModel.remark_stack_outperform_svd_stack
#print axioms StackedSVD.MultiTableModel.remark_stack_outperform_svd_svdstack
#print axioms StackedSVD.MultiTableModel.remark_svd_outperform_stack_pair_stack
#print axioms StackedSVD.MultiTableModel.remark_svd_outperform_stack_pair_svdstack
#print axioms StackedSVD.MultiTableModel.remark_svd_outperform_stack_three_svdstack
#print axioms StackedSVD.MultiTableModel.remark_svd_outperform_stack_three_stack
#print axioms StackedSVD.MultiTableModel.remark_svd_outperform_stack_three_binary
#print axioms StackedSVD.MultiTableModel.remark_svd_outperform_stack_two

/-! ### Uniform Rayleigh bound (P1) -/

#print axioms StackedSVD.MultiTableModel.svdstackPerfW_le_rowBound
#print axioms StackedSVD.MultiTableModel.svdstackPerfW_uniform_bound
#print axioms StackedSVD.MultiTableModel.thm_svdstack_weighted_gaussian_opt_full

/-! ### `app:wstacksvd_mle` -/

#print axioms StackedSVD.MultiTableModel.mleLogLik_eq
#print axioms StackedSVD.MultiTableModel.mleLogLik_le_of_mem_topSpace
#print axioms StackedSVD.MultiTableModel.mleLogLik_le_iff

/-! ### `lem:secular_equation` -/

#print axioms StackedSVD.Secular.card_sub_one_le_finrank_eigenspace
#print axioms StackedSVD.Secular.secularDiag_stack

/-! ## Part 3: rank r (added 2026-09-02 for AUDIT_PACK_V4; Track C stage 1 only) -/

/-! ### Part 1: rank-1 names missing from `scripts/AuditSignatures.lean` -/

#check @StackedSVD.MultiTableModel.rowBound_tendsto

/-! ### Part 2: rank r, models and hypothesis structures -/

#check @StackedSVD.UnalignedModel
#check @StackedSVD.UnalignedModel.SubspaceLaw
#check @StackedSVD.SpikedModelR
#check @StackedSVD.SpikedModelR.TableLawR
#check @StackedSVD.UnalignedModelR
#check @StackedSVD.UnalignedModelR.SubspaceLawG
#check @StackedSVD.UnalignedModelR.IndepNoise
#check @StackedSVD.RankRStack

/-! ### Part 3: Track A, rank-r stacksvd for Gaussian noise -/

#check @StackedSVD.UnalignedModel.prop_stacksvd_subspace
#check @StackedSVD.RankRStack.align_of_gaussian
#check @StackedSVD.UnalignedModel.subspaceLaw_of_gaussian
#check @StackedSVD.UnalignedModel.prop_stacksvd_subspace_gaussian
#check @StackedSVD.UnalignedModel.topGap_stackGram_whp_of_gaussian
#check @StackedSVD.UnalignedModel.prop_stacksvd_subspace_frobenius_eig_gaussian
#check @StackedSVD.UnalignedModelR.subspaceLawG_of_gaussian
#check @StackedSVD.UnalignedModelR.prop_stacksvd_subspace_general_gaussian
#check @StackedSVD.RankR.Example.example_perfStackR_tendsto_gaussian

/-! ### Part 4: Track B, svdstack at general `r_i` -/

#check @StackedSVD.UnalignedModelR.lem_general_rank_delocalization_general
#check @StackedSVD.UnalignedModelR.gramR_general
#check @StackedSVD.UnalignedModelR.VtV_tendsto_general
#check @StackedSVD.UnalignedModelR.prop_general_rank_unweighted_svdstack_general
#check @StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_conv
#check @StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r
#check @StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_paper
#check @StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_full
#check @StackedSVD.UnalignedModelR.perfRGW_uniform_bound
#check @StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_full_gaussian_one
#check @StackedSVD.UnalignedModelR.prop_general_rank_unweighted_svdstack_general_gaussian_one
#check @StackedSVD.UnalignedModelR.prop_general_rank_unweighted_svdstack_general_frobenius_eig
#check @StackedSVD.SpikedModelR.tableLawR_of_gaussian
#check @StackedSVD.UnalignedModelR.tableLawR_of_gaussian
#check @StackedSVD.UnalignedModel.thm_gen_rank_weight_svdstak_full
#check @StackedSVD.UnalignedModel.thm_gen_rank_weight_svdstak_full_gaussian
#check @StackedSVD.UnalignedModel.perfRW_uniform_bound
#check @StackedSVD.UnalignedModel.prop_general_rank_unweighted_svdstack_gaussian
#check @StackedSVD.UnalignedModel.thm_gen_rank_weight_svdstak_gaussian

/-! ### Part 5: Track D item D3, the aggregate clause -/

#check @StackedSVD.Sagg
#check @StackedSVD.trace_conj_inv_ABlock_aligned
#check @StackedSVD.limitOptG_aligned
#check @StackedSVD.UnalignedModelR.thm_rank_r_svdstack_aggregate

/-! ### Part 6: Track C stage 1 -/

#check @StackedSVD.ScalarsC.rhoSq_lt_rhoSq
#check @StackedSVD.RankRStack.tendsto_measure_lamMax_w1R_le_tau
#check @StackedSVD.RankRStack.tendsto_measure_eigenvalues₀_le_tau
#check @StackedSVD.SpikedModelR.toStack
#check @StackedSVD.SpikedModelR.exists_perm_coreEig
#check @StackedSVD.Frame.lt_eigenvalues₀_iff_of_card
#check @StackedSVD.specProjIdx_eq_specProj_Ioc
#check @StackedSVD.overlapIdx_eq_normSq_specProj_Ioi_sub

/-! ### Part 7: axiom sets -/

#print axioms StackedSVD.MultiTableModel.svdstackPerfW_le_opt_whp
#print axioms StackedSVD.MultiTableModel.rowBound_tendsto
#print axioms StackedSVD.svdstackLimitOpt_eq_of_single
#print axioms StackedSVD.svdstackLimit_pair
#print axioms StackedSVD.remark_three_stack_tendsto_zero
#print axioms StackedSVD.svdstackLimitOpt_zero
#print axioms StackedSVD.MultiTableModel.remark_stack_outperform_svd_svdstack_uniform
#print axioms StackedSVD.MultiTableModel.remark_svd_outperform_stack_two_svdstack_unweighted
#print axioms StackedSVD.inner_toOp_self_eq_lamMax_iff
#print axioms StackedSVD.MultiTableModel.mem_topSpace_of_mleLogLik_max
#print axioms StackedSVD.MultiTableModel.mleLogLik_max_iff_mem_topSpace

#print axioms StackedSVD.UnalignedModel.prop_stacksvd_subspace
#print axioms StackedSVD.RankRStack.align_of_gaussian
#print axioms StackedSVD.UnalignedModel.subspaceLaw_of_gaussian
#print axioms StackedSVD.UnalignedModel.prop_stacksvd_subspace_gaussian
#print axioms StackedSVD.UnalignedModel.topGap_stackGram_whp_of_gaussian
#print axioms StackedSVD.UnalignedModel.prop_stacksvd_subspace_frobenius_eig_gaussian
#print axioms StackedSVD.UnalignedModelR.subspaceLawG_of_gaussian
#print axioms StackedSVD.UnalignedModelR.prop_stacksvd_subspace_general_gaussian
#print axioms StackedSVD.RankR.Example.example_perfStackR_tendsto_gaussian

#print axioms StackedSVD.UnalignedModelR.lem_general_rank_delocalization_general
#print axioms StackedSVD.UnalignedModelR.gramR_general
#print axioms StackedSVD.UnalignedModelR.VtV_tendsto_general
#print axioms StackedSVD.UnalignedModelR.prop_general_rank_unweighted_svdstack_general
#print axioms StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_conv
#print axioms StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r
#print axioms StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_paper
#print axioms StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_full
#print axioms StackedSVD.UnalignedModelR.perfRGW_uniform_bound
#print axioms StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_full_gaussian_one
#print axioms StackedSVD.UnalignedModelR.prop_general_rank_unweighted_svdstack_general_gaussian_one
#print axioms StackedSVD.UnalignedModelR.prop_general_rank_unweighted_svdstack_general_frobenius_eig
#print axioms StackedSVD.SpikedModelR.tableLawR_of_gaussian
#print axioms StackedSVD.UnalignedModelR.tableLawR_of_gaussian
#print axioms StackedSVD.UnalignedModel.thm_gen_rank_weight_svdstak_full
#print axioms StackedSVD.UnalignedModel.thm_gen_rank_weight_svdstak_full_gaussian
#print axioms StackedSVD.UnalignedModel.perfRW_uniform_bound
#print axioms StackedSVD.UnalignedModel.prop_general_rank_unweighted_svdstack_gaussian
#print axioms StackedSVD.UnalignedModel.thm_gen_rank_weight_svdstak_gaussian

#print axioms StackedSVD.trace_conj_inv_ABlock_aligned
#print axioms StackedSVD.limitOptG_aligned
#print axioms StackedSVD.UnalignedModelR.thm_rank_r_svdstack_aggregate

#print axioms StackedSVD.ScalarsC.rhoSq_lt_rhoSq
#print axioms StackedSVD.RankRStack.tendsto_measure_lamMax_w1R_le_tau
#print axioms StackedSVD.RankRStack.tendsto_measure_eigenvalues₀_le_tau
#print axioms StackedSVD.SpikedModelR.exists_perm_coreEig
#print axioms StackedSVD.Frame.lt_eigenvalues₀_iff_of_card
#print axioms StackedSVD.specProjIdx_eq_specProj_Ioc
#print axioms StackedSVD.overlapIdx_eq_normSq_specProj_Ioi_sub

/-! ## Part 4: rank r, Track C stages 2 to 4 and Track D (added 2026-09-02 for
AUDIT_PACK_V4.1). Same convention as Part 3: signatures first, one axiom-set subsection last. -/

/-! ### Part 4.1: Track C stage 2 (task C3 `RankR/RMT/AlignTauR.lean`, task C6
`RankR/RMT/DelocAffineR.lean`, task C9 `RankR/RMT/R6R.lean`) -/

#check @StackedSVD.OutliersR.supercritical_of_lt_rhoSq
#check @StackedSVD.OutliersR.exists_sum_equiv_split_rho
#check @StackedSVD.RankRStack.tendstoInProb_normSq_specProj_Ioi_tau
#check @StackedSVD.RankRStack.tendsto_measure_count_Ioi_tau

#check @StackedSVD.overlapIdx_mul_right
#check @StackedSVD.mul_right_of_transpose_mulVec_col
#check @StackedSVD.sum_overlapIdx_le_one
#check @StackedSVD.lintegral_overlapIdx_eq_of_orth
#check @StackedSVD.lintegral_overlapIdx_affine_le
#check @StackedSVD.SpikedModelR.lintegral_overlapIdx_model_le
#check @StackedSVD.SpikedModelR.measure_overlapIdx_ge_le
#check @StackedSVD.SpikedModelR.delocUniform_of_gaussian

#check @StackedSVD.normSq_specProjIdx_le_edge
#check @StackedSVD.SpikedModelR.exists_split_subcritical
#check @StackedSVD.SpikedModelR.align_cross_of_gaussian_subcritical_aux
#check @StackedSVD.SpikedModelR.align_cross_of_gaussian_subcritical

/-! ### Part 4.2: Track C stage 3 (task C7 `RankR/RMT/TableSimple.lean`, task C5
`RankR/RMT/TableAlignSup.lean`) -/

#check @StackedSVD.SpikedModelR.simple_of_gaussian

#check @StackedSVD.RankRStack.measurableSet_count_Ioi
#check @StackedSVD.SpikedModelR.rhoSq_lt_rhoSq_of_lt
#check @StackedSVD.SpikedModelR.rhoSq_lt_rhoSq_of_gt
#check @StackedSVD.SpikedModelR.exists_margin
#check @StackedSVD.SpikedModelR.card_filter_rhoSq_coreEig
#check @StackedSVD.SpikedModelR.eigenvalues₀_toStack_gram
#check @StackedSVD.SpikedModelR.align_cross_of_gaussian_supercritical_aux
#check @StackedSVD.SpikedModelR.align_cross_of_gaussian_supercritical

/-! ### Part 4.3: Track C stage 4 (task C10, `RankR/RMT/TableLawGaussian.lean`, the
unconditional Gaussian `TableLawR`) and its Section 7 consumers in `RankR/GeneralGaussian.lean`
(the seven `_gaussian` twins; names only, from `notes/INTERFACES.md`, no source read of that
file per the coordinator's concurrency rule) -/

#check @StackedSVD.SpikedModelR.tableLawR_of_gaussian_rk
#check @StackedSVD.UnalignedModelR.tableLawR_of_gaussian_rk

#check @StackedSVD.UnalignedModelR.prop_general_rank_unweighted_svdstack_general_gaussian
#check @StackedSVD.UnalignedModelR.prop_general_rank_unweighted_svdstack_general_frobenius_eig_gaussian
#check @StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_frobenius_eig_gaussian
#check @StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_of_rank_gaussian
#check @StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_paper_gaussian
#check @StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_full_gaussian
#check @StackedSVD.UnalignedModelR.thm_rank_r_svdstack_aggregate_gaussian

/-! ### Part 4.4: Track D item D1 (`RankR/StackGamma.lean`, the definitions; the hypothesis
structure `HeteroLawR` and `thm_rank_r_stacksvd` were not in the tree when this part was
written, see `notes/archive/rankr_D1_statement.md`. They landed with Track D and Track E; Part 5
below carries them and their Gaussian discharge.) -/

#check @StackedSVD.Scalars.wStackR
#check @StackedSVD.Scalars.thetaTildeSq
#check @StackedSVD.Scalars.gammaR
#check @StackedSVD.Scalars.ellR
#check @StackedSVD.UnalignedModelR.stackXW
#check @StackedSVD.UnalignedModelR.thetaAligned
#check @StackedSVD.UnalignedModelR.ellR_thetaAligned
#check @StackedSVD.UnalignedModelR.stackXJ
#check @StackedSVD.UnalignedModelR.stackGramJ
#check @StackedSVD.UnalignedModelR.stackOverlapJ
#check @StackedSVD.UnalignedModelR.vhatStackR
#check @StackedSVD.UnalignedModelR.frobSqStackR

/-! ### Part 4.5: Track D item D4-lin (`LinAlg/SpecIdxPerturb.lean`, index-`j` projector
continuity) -/

#check @StackedSVD.eigVal
#check @StackedSVD.continuousAt_eigVal_symMat
#check @StackedSVD.normSq_specProjIdx_eq_sub
#check @StackedSVD.normSq_specProjTop_continuousAt
#check @StackedSVD.compFun
#check @StackedSVD.compFun_eq
#check @StackedSVD.continuousAt_compFun
#check @StackedSVD.SimpleIdx
#check @StackedSVD.normSq_specProjIdx_eq_inner_sq
#check @StackedSVD.specIdxSimple_whp_of_tendsto

/-! ### Part 4.6: Track D item D4-spec (`RankR/AlignedComponent.lean`, the aligned limit
spectrum at index `j`, and the three lemmas `LinAlg/Eigen.lean` gained for it) -/

#check @StackedSVD.BBlockWcol
#check @StackedSVD.ABlockW_optWG_aligned_eq
#check @StackedSVD.transpose_BBlockW_optWG_aligned
#check @StackedSVD.eigenvalues₀_ABlockW_optWG_aligned
#check @StackedSVD.SaggSep
#check @StackedSVD.Sagg_sep_of_strictAnti
#check @StackedSVD.topGap_ABlockW_optWG_aligned_of_sep
#check @StackedSVD.topGap_succ_ABlockW_optWG_aligned_of_sep
#check @StackedSVD.simpleSpec_ABlockW_optWG_aligned_of_sep
#check @StackedSVD.normSq_specProjIdx_ABlockW_aligned_of_sep
#check @StackedSVD.compLimit_aligned_of_sep
#check @StackedSVD.compLimit_aligned
#check @StackedSVD.Sagg_antitone
#check @StackedSVD.Sagg_lt_of_pos
#check @StackedSVD.UnalignedModelR.saggSep_of_model

#check @StackedSVD.antitone_padSpec
#check @StackedSVD.prod_X_sub_C_padSpec
#check @StackedSVD.eigenvalues₀_mul_transpose_of_transpose_mul_diagonal

/-! ### Part 4.7: axiom sets for Part 4.1 to 4.6 -/

#print axioms StackedSVD.OutliersR.supercritical_of_lt_rhoSq
#print axioms StackedSVD.OutliersR.exists_sum_equiv_split_rho
#print axioms StackedSVD.RankRStack.tendstoInProb_normSq_specProj_Ioi_tau
#print axioms StackedSVD.RankRStack.tendsto_measure_count_Ioi_tau

#print axioms StackedSVD.overlapIdx_mul_right
#print axioms StackedSVD.mul_right_of_transpose_mulVec_col
#print axioms StackedSVD.sum_overlapIdx_le_one
#print axioms StackedSVD.lintegral_overlapIdx_eq_of_orth
#print axioms StackedSVD.lintegral_overlapIdx_affine_le
#print axioms StackedSVD.SpikedModelR.lintegral_overlapIdx_model_le
#print axioms StackedSVD.SpikedModelR.measure_overlapIdx_ge_le
#print axioms StackedSVD.SpikedModelR.delocUniform_of_gaussian

#print axioms StackedSVD.normSq_specProjIdx_le_edge
#print axioms StackedSVD.SpikedModelR.exists_split_subcritical
#print axioms StackedSVD.SpikedModelR.align_cross_of_gaussian_subcritical_aux
#print axioms StackedSVD.SpikedModelR.align_cross_of_gaussian_subcritical

#print axioms StackedSVD.SpikedModelR.simple_of_gaussian

#print axioms StackedSVD.RankRStack.measurableSet_count_Ioi
#print axioms StackedSVD.SpikedModelR.rhoSq_lt_rhoSq_of_lt
#print axioms StackedSVD.SpikedModelR.rhoSq_lt_rhoSq_of_gt
#print axioms StackedSVD.SpikedModelR.exists_margin
#print axioms StackedSVD.SpikedModelR.card_filter_rhoSq_coreEig
#print axioms StackedSVD.SpikedModelR.eigenvalues₀_toStack_gram
#print axioms StackedSVD.SpikedModelR.align_cross_of_gaussian_supercritical_aux
#print axioms StackedSVD.SpikedModelR.align_cross_of_gaussian_supercritical

#print axioms StackedSVD.SpikedModelR.tableLawR_of_gaussian_rk
#print axioms StackedSVD.UnalignedModelR.tableLawR_of_gaussian_rk

#print axioms StackedSVD.UnalignedModelR.prop_general_rank_unweighted_svdstack_general_gaussian
#print axioms StackedSVD.UnalignedModelR.prop_general_rank_unweighted_svdstack_general_frobenius_eig_gaussian
#print axioms StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_frobenius_eig_gaussian
#print axioms StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_of_rank_gaussian
#print axioms StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_paper_gaussian
#print axioms StackedSVD.UnalignedModelR.thm_gen_rank_weight_svdstak_general_r_full_gaussian
#print axioms StackedSVD.UnalignedModelR.thm_rank_r_svdstack_aggregate_gaussian

#print axioms StackedSVD.Scalars.wStackR
#print axioms StackedSVD.Scalars.thetaTildeSq
#print axioms StackedSVD.Scalars.gammaR
#print axioms StackedSVD.Scalars.ellR
#print axioms StackedSVD.UnalignedModelR.stackXW
#print axioms StackedSVD.UnalignedModelR.thetaAligned
#print axioms StackedSVD.UnalignedModelR.ellR_thetaAligned
#print axioms StackedSVD.UnalignedModelR.stackXJ
#print axioms StackedSVD.UnalignedModelR.stackGramJ
#print axioms StackedSVD.UnalignedModelR.stackOverlapJ
#print axioms StackedSVD.UnalignedModelR.vhatStackR
#print axioms StackedSVD.UnalignedModelR.frobSqStackR

#print axioms StackedSVD.eigVal
#print axioms StackedSVD.continuousAt_eigVal_symMat
#print axioms StackedSVD.normSq_specProjIdx_eq_sub
#print axioms StackedSVD.normSq_specProjTop_continuousAt
#print axioms StackedSVD.compFun
#print axioms StackedSVD.compFun_eq
#print axioms StackedSVD.continuousAt_compFun
#print axioms StackedSVD.SimpleIdx
#print axioms StackedSVD.normSq_specProjIdx_eq_inner_sq
#print axioms StackedSVD.specIdxSimple_whp_of_tendsto

#print axioms StackedSVD.BBlockWcol
#print axioms StackedSVD.ABlockW_optWG_aligned_eq
#print axioms StackedSVD.transpose_BBlockW_optWG_aligned
#print axioms StackedSVD.eigenvalues₀_ABlockW_optWG_aligned
#print axioms StackedSVD.SaggSep
#print axioms StackedSVD.Sagg_sep_of_strictAnti
#print axioms StackedSVD.topGap_ABlockW_optWG_aligned_of_sep
#print axioms StackedSVD.topGap_succ_ABlockW_optWG_aligned_of_sep
#print axioms StackedSVD.simpleSpec_ABlockW_optWG_aligned_of_sep
#print axioms StackedSVD.normSq_specProjIdx_ABlockW_aligned_of_sep
#print axioms StackedSVD.compLimit_aligned_of_sep
#print axioms StackedSVD.compLimit_aligned
#print axioms StackedSVD.Sagg_antitone
#print axioms StackedSVD.Sagg_lt_of_pos
#print axioms StackedSVD.UnalignedModelR.saggSep_of_model

#print axioms StackedSVD.antitone_padSpec
#print axioms StackedSVD.prod_X_sub_C_padSpec
#print axioms StackedSVD.eigenvalues₀_mul_transpose_of_transpose_mul_diagonal

/-! ## Part 5: Track E, the rank-`r` weighted stackSVD endpoint (added 2026-09-02)

`thm:rank_r_stacksvd` (`main_paper.tex:2337`, Appendix E, `sec:rank_r`). Layer 1 takes the
hypothesis structure `HeteroLawR`; the Gaussian facade discharges it, so
`thm_rank_r_stacksvd_gaussian` has model hypotheses only. The model is
`UnalignedModelR μ M n d r (alignedRk M r)`, the exactly aligned class, with `hθpos` and
`hθanti` on each table (`notes/FLAGGED.md` item 17 (7), `notes/paper_edits.md` finding E7). -/

/-! ### Part 5.1: the hypothesis structures -/

#check @StackedSVD.UnalignedModelR.HeteroLawR
#check @StackedSVD.UnalignedModelR.HeteroEdgeR

/-! ### Part 5.2: Layer 1, the four clauses of `thm:rank_r_stacksvd` -/

#check @StackedSVD.UnalignedModelR.thm_rank_r_stacksvd
#check @StackedSVD.UnalignedModelR.thm_rank_r_stacksvd_proj
#check @StackedSVD.UnalignedModelR.thm_rank_r_stacksvd_inner
#check @StackedSVD.UnalignedModelR.thm_rank_r_stacksvd_frobenius

/-! ### Part 5.3: the Gaussian discharge -/

#check @StackedSVD.UnalignedModelR.heteroEdgeR_of_gaussian
#check @StackedSVD.UnalignedModelR.heteroLawR_of_gaussian
#check @StackedSVD.UnalignedModelR.thm_rank_r_stacksvd_gaussian
#check @StackedSVD.UnalignedModelR.thm_rank_r_stacksvd_proj_gaussian
#check @StackedSVD.UnalignedModelR.thm_rank_r_stacksvd_inner_gaussian
#check @StackedSVD.UnalignedModelR.thm_rank_r_stacksvd_frobenius_gaussian

/-! ### Part 5.4: axiom sets for Part 5 -/

#print axioms StackedSVD.UnalignedModelR.thm_rank_r_stacksvd
#print axioms StackedSVD.UnalignedModelR.heteroEdgeR_of_gaussian
#print axioms StackedSVD.UnalignedModelR.heteroLawR_of_gaussian
#print axioms StackedSVD.UnalignedModelR.thm_rank_r_stacksvd_gaussian

/-! ## Part 6: the satisfiability witnesses (added 2026-09-02)

`StackedSVD/Sat.lean` builds concrete Gaussian models and applies the three headline Gaussian
endpoints to them, with no hypothesis left over. The three theorems below are therefore the
evidence that the hypothesis sets are satisfiable and the endpoints are not vacuous. -/

#check @StackedSVD.Sat.Rank1.sat_stacksvd_weighted
#check @StackedSVD.Sat.Rank1.sat_svdstack_weighted
#check @StackedSVD.Sat.RankR.sat_rank_r_stacksvd

#print axioms StackedSVD.Sat.Rank1.sat_stacksvd_weighted
#print axioms StackedSVD.Sat.Rank1.sat_svdstack_weighted
#print axioms StackedSVD.Sat.RankR.sat_rank_r_stacksvd
