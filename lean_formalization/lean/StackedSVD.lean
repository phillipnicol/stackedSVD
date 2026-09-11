import StackedSVD.Defs
import StackedSVD.RMT
import StackedSVD.Spectral
import StackedSVD.LinAlg.TopProjPerturb
import StackedSVD.LinAlg.SpecProjPerturb
import StackedSVD.LinAlg.KyFan
import StackedSVD.LinAlg.Eigen
import StackedSVD.LinAlg.Frame
import StackedSVD.Prob.TendstoInProb
import StackedSVD.Prob.PolynomialNull
import StackedSVD.Prob.GaussianMatrix
import StackedSVD.Prob.GaussianAdapters
import StackedSVD.Prob.WithDensityPi
import StackedSVD.Prob.GaussianDensity
import StackedSVD.Prob.NoiseLaw
import StackedSVD.Prob.LinFormMoments
import StackedSVD.Prob.NoiseMoments
-- Stage 1 of the non-Gaussian extension (notes/archive/prop_single_table_general.md): the two
-- shared probability lemmas of units G1b and G0.
import StackedSVD.Prob.Tensorization
import StackedSVD.Prob.Chebyshev
import StackedSVD.StackSVD
import StackedSVD.StackSVD.Main
import StackedSVD.Scalars
import StackedSVD.StackSVDWeighted
import StackedSVD.StackSVD.Weighted
import StackedSVD.Secular
import StackedSVD.ThetaEst
import StackedSVD.SVDStack
import StackedSVD.SVDStack.DelocDir
import StackedSVD.SVDStack.DelocFinite
import StackedSVD.StrictFacades
import StackedSVD.MLE
import StackedSVD.Remarks
import StackedSVD.MLEConverse
import StackedSVD.MLEMarginal.Defs
import StackedSVD.MLEMarginal.RowLaw
import StackedSVD.MLEMarginal.TableLaw
import StackedSVD.MLEMarginal.LogLik
import StackedSVD.MLEMarginal.Main
import StackedSVD.RemarksUniform
import StackedSVD.RankR.Defs
import StackedSVD.RankR.Unweighted
import StackedSVD.RankR.Weighted
import StackedSVD.RankR.WeightedMain
import StackedSVD.RankR.Example
import StackedSVD.RankR.Subspace
import StackedSVD.RankR.Frobenius
import StackedSVD.RankR.General
import StackedSVD.RankR.Flatten
-- Layer 2 (Gaussian proof of SingleTableLaw), by roadmap item
import StackedSVD.RMT.MP
import StackedSVD.RMT.MP7
import StackedSVD.RMT.R0
import StackedSVD.RMT.R3
import StackedSVD.RMT.R3minus
import StackedSVD.RMT.R6
import StackedSVD.RMT.T
import StackedSVD.RMT.ResolvDeriv
import StackedSVD.RMT.R1
import StackedSVD.RMT.R2
import StackedSVD.RMT.R5
import StackedSVD.RMT.Sup
import StackedSVD.RMT.TailShift
import StackedSVD.RMT.Full
import StackedSVD.RMT.R4
import StackedSVD.RMT.R4C
import StackedSVD.RMT.Simplicity
import StackedSVD.RMT.Symmetry
import StackedSVD.RMT.Het.MPhet
-- Layer 2 heteroscedastic (Gaussian proof of HeteroLaw, notes/archive/plan_heterolaw_A.md)
import StackedSVD.RMT.Het.Split
import StackedSVD.RMT.Het.Duality
import StackedSVD.RMT.Het.R4het
import StackedSVD.RMT.Het.R3het
import StackedSVD.RMT.Het.Stein
import StackedSVD.RMT.Het.R1het
import StackedSVD.RMT.Het.R2het
import StackedSVD.RMT.Het.R5het
import StackedSVD.RMT.Het.Simplicity
import StackedSVD.RMT.Het.R6het
import StackedSVD.RMT.Het.EdgeScalar
import StackedSVD.RMT.Het.EdgeSharp
import StackedSVD.RMT.Het.Sup
-- Layer 2 rank-r (Gaussian proof of SubspaceLaw, notes/archive/rankr_plan_A.md)
import StackedSVD.RankR.RMT.Stack
import StackedSVD.RankR.RMT.Split
import StackedSVD.RankR.RMT.Forms
import StackedSVD.RankR.RMT.Outliers
import StackedSVD.RankR.RMT.SimplicityR
import StackedSVD.RankR.RMT.DelocR
import StackedSVD.RankR.RMT.EdgeR
import StackedSVD.RankR.RMT.SymmetryR
import StackedSVD.RankR.RMT.DecompR
import StackedSVD.RankR.RMT.SimplicityAffineR
import StackedSVD.RankR.RMT.EdgeDetR
import StackedSVD.RankR.RMT.OutliersG
import StackedSVD.RankR.RMT.AlignOutG
import StackedSVD.RankR.RMT.ShiftR
import StackedSVD.RankR.RMT.EdgeGlueDetR
import StackedSVD.RankR.RMT.EdgeGlueR
import StackedSVD.RankR.RMT.AlignG
import StackedSVD.RankR.RMT.EdgeTauR
import StackedSVD.RankR.RMT.TableStack
import StackedSVD.RankR.RMT.DelocAffineR
import StackedSVD.RankR.RMT.R6R
import StackedSVD.RankR.RMT.AlignTauR
import StackedSVD.RankR.RMT.TableSimple
import StackedSVD.RankR.RMT.TableAlignSup
import StackedSVD.RankR.RMT.TableLawGaussian
import StackedSVD.LinAlg.SpecIdxMeas
import StackedSVD.LinAlg.SpecWindow
import StackedSVD.LinAlg.SpecIdxPerturb
import StackedSVD.RankR.GramR
import StackedSVD.RankR.GeneralMain
import StackedSVD.RankR.GeneralFrob
import StackedSVD.RankR.SubspaceG
import StackedSVD.RankR.SubspaceMain
import StackedSVD.RankR.SubspaceGStack
import StackedSVD.RankR.SubspaceGaussian
import StackedSVD.RankR.Aligned
import StackedSVD.RankR.AlignedComponent
import StackedSVD.RankR.AlignedMain
import StackedSVD.RankR.WeightedUpper
import StackedSVD.RankR.WeightedUpperG
import StackedSVD.RankR.GeneralGaussian
import StackedSVD.RankR.AlignedOrth
import StackedSVD.RankR.SingleWeight.Scalars
import StackedSVD.RankR.SingleWeight.ScalarsOne
import StackedSVD.RankR.SingleWeight.Defs
import StackedSVD.RankR.SingleWeight.Main
import StackedSVD.RankR.SingleWeight.Example
import StackedSVD.RankR.SingleWeight.Existence
import StackedSVD.RankR.SingleWeight.Het.Scalars
import StackedSVD.RankR.SingleWeight.Het.Forms
import StackedSVD.RankR.SingleWeight.Het.Frame
import StackedSVD.RankR.SingleWeight.Het.Count
import StackedSVD.RankR.SingleWeight.Het.Outliers
import StackedSVD.RankR.SingleWeight.Het.Align
import StackedSVD.RankR.SingleWeight.Het.Sub
import StackedSVD.RankR.SingleWeight.Het.Sup
import StackedSVD.RankR.SingleWeight.Regimes
import StackedSVD.RankR.SingleWeight.Tie
import StackedSVD.RankR.SingleWeight.Het.OneOutDet
import StackedSVD.RankR.SingleWeight.Het.OneOutCount
import StackedSVD.RankR.SingleWeight.Het.OneOutAlign
import StackedSVD.RankR.SingleWeight.Het.OneOutBulk
import StackedSVD.RankR.SingleWeight.Suboptimality
import StackedSVD.RankR.Het.Scalars
import StackedSVD.RankR.Het.Split
import StackedSVD.RankR.Het.Duality
import StackedSVD.RankR.Het.Simplicity
import StackedSVD.RankR.Het.Edge
import StackedSVD.RankR.Het.Forms
import StackedSVD.RankR.Het.Deloc
import StackedSVD.RankR.Het.BulkDet
import StackedSVD.RankR.Het.Outliers
import StackedSVD.RankR.Het.Align
import StackedSVD.RankR.Het.Bulk
import StackedSVD.RankR.Het.Sup
import StackedSVD.RankR.StackGamma
import StackedSVD.RankR.StackMain
-- The satisfiability witnesses: concrete Gaussian models that meet every hypothesis of the
-- three headline Gaussian endpoints. Imported so the axiom gate audits them.
import StackedSVD.Sat
import StackedSVD.Existence
-- Vendored COLT83 (RemyDegenne/colt-2026-83, Apache-2.0): Stein identity, Gaussian
-- interpolation, Sudakov-Fernique, Borell-TIS. See Vendor/COLT83/README.md.
import StackedSVD.Vendor.COLT83.Mathlib.Probability.SteinReal
import StackedSVD.Vendor.COLT83.Mathlib.Probability.SteinIdentity
import StackedSVD.Vendor.COLT83.Mathlib.Probability.GaussianInterpolation
import StackedSVD.Vendor.COLT83.Mathlib.Probability.SudakovFernique
import StackedSVD.Vendor.COLT83.Mathlib.Probability.BorellTIS
-- Stage 1 of the non-Gaussian extension (notes/archive/prop_single_table_general.md): the
-- proof units of the general-law single table, the structure `SpikedModel.ResolventFormsC`
-- (`RMT/General/Defs.lean`; moved there from `RMT/General/Statements.lean` by follow-up item
-- F35, 2026-09-09, which retired that file, the file that had held the 12 target statements
-- of the plan while they carried `sorry`), the endpoint (`Sup.lean`) and the Layer 1 facades
-- (`General/Layer1.lean`). The group is ordered by dependency, deepest first.
import StackedSVD.RMT.General.Defs
import StackedSVD.RMT.General.QuadForm
import StackedSVD.RMT.General.Stability
import StackedSVD.RMT.General.Companion
import StackedSVD.RMT.General.MPtilde
import StackedSVD.RMT.General.Trace
import StackedSVD.RMT.General.R3minus
import StackedSVD.RMT.General.Simplicity
import StackedSVD.RMT.General.Iso
import StackedSVD.RMT.General.IsoMixed
import StackedSVD.RMT.General.FormsBridge
import StackedSVD.RMT.General.ProbC
import StackedSVD.RMT.General.Forms
import StackedSVD.RMT.General.FormsLimits
import StackedSVD.RMT.General.Deloc
import StackedSVD.RMT.General.FormsSup
import StackedSVD.RMT.General.FormsGeneral
import StackedSVD.RMT.General.DelocLimits
import StackedSVD.RMT.General.DelocAlign
import StackedSVD.RMT.General.DelocUniform
import StackedSVD.RMT.General.DelocUniformSup
import StackedSVD.RMT.General.Sup
-- Stage 3 of the non-Gaussian extension (notes/stage3_edge.md, notes/STAGE3_CAMPAIGN.md): the
-- sharp upper edge at a general law (the upper half of Bai-Yin at four moments), route C, the
-- moment method with a truncation. `Edge/Sup.lean` holds the endpoint
-- `SpikedModel.opNorm_sq_edge_of_general` and the corollary `lamMax_W0_edge_of_general`, the
-- `hedge` hypothesis of `singleTableLaw_of_general` verbatim. Deepest first.
import StackedSVD.RMT.General.Edge.Defs
import StackedSVD.RMT.General.Edge.Arith
import StackedSVD.RMT.General.Edge.Trunc
import StackedSVD.RMT.General.Edge.Sparse
import StackedSVD.RMT.General.Edge.Trace
import StackedSVD.RMT.General.Edge.Compare
import StackedSVD.RMT.General.Edge.Gaussian
import StackedSVD.RMT.General.Edge.Count
import StackedSVD.RMT.General.Edge.Code
import StackedSVD.RMT.General.Edge.Dyck
import StackedSVD.RMT.General.Edge.CountBound
import StackedSVD.RMT.General.Edge.Excess
import StackedSVD.RMT.General.Edge.Markov
import StackedSVD.RMT.General.Edge.Sup
import StackedSVD.General.Layer1
-- The statement file: every result of the paper, in the order of docs/THEOREMS.md. Placed
-- last because it imports across nearly the whole tree (this file is not alphabetically
-- ordered; it groups files by layer and dependency, and Main.lean depends on the layers
-- above).
import StackedSVD.Main
