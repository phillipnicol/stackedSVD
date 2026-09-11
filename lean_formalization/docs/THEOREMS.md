# THEOREMS.md: what is proved, and under which hypotheses

> Location (2026-09-03): this file lives in `docs/`. The bare names `TECHNICAL.md`, `THEOREMS.md`
> and `AUDIT_DOC.md` refer to the files of this folder; every other path is relative to the
> repository root. `paper_edits.md` and `NOTE_FOR_COAUTHORS.md` moved from this folder to
> `notes/` on 2026-09-11; this file names them `notes/paper_edits.md` and
> `notes/NOTE_FOR_COAUTHORS.md` throughout.
> The short `README.md` at the root is the plain-language entry point.

This document is self-contained. It states every result of the paper *"Stacked SVD or SVD
stacked? A random matrix theory perspective on data integration"* that the Lean development
covers, next to the Lean declaration that carries it. An auditor can read it without the
repository open.

Source of the paper text: `main_paper.tex`, the read-only snapshot in `notes/paper/` of the
development repository (taken 2026-08-29 from the AoS resubmission; not part of this public
copy, and the arXiv version carries the same labels). Label references are labels, not line numbers;
grep the label.

Source of the Lean text: `lean/StackedSVD/`, branch `stage3`. Section 8 names the
commit of the build of record. Every signature below was copied from the tree (the first
rows on 2026-09-02, later rows with their own dates) and is re-checked against the build of
record by `python3 scripts/check_theorems_sigs.py docs/THEOREMS.md` (gate 4; section 8 holds
the last run).
The declaration names never change; a redundant hypothesis may be dropped from a signature by
a later cleanup pass, and `notes/archive/cleanup_2026-09-02_removed.md` records every such
change.

## 0. How to read a result section

Each section has up to eight parts.

1. **Paper.** The theorem environment, verbatim LaTeX.
2. **Paper, as amended.** A result that a paper edit touches carries this second block: the
   whole statement environment, verbatim as it reads after the replacement, or a one-line note
   when the edit changes only prose outside the environment. The exact replacements, old text
   and new text, are the **Exact replacement** subsections of `notes/paper_edits.md`.
   `python3 scripts/check_paper_edits.py` checks that every old text occurs exactly once in the
   snapshot `main_paper.tex` and that every amended block here is verbatim in the amended paper.
3. **Lean.** The `theorem` signature, verbatim, for the Gaussian form and for the Layer 1
   form when both exist, and, in the five sections that have one, a **Lean, general noise**
   block for the four-moment form (3.1, 3.2, 3.4, 3.5 and 5.4); section 2.1 quotes the same
   general-noise discharge next to the Gaussian one. Exceptions, for readability: a
   proof-term argument (a `by simpa
   ...` witness) is written `…` or `_`; an implicit binder or a binder type that the
   surrounding text explains may be dropped; an `eigenvalues₀` index may be paraphrased; a
   docstring may be shortened; a definition block may collect definitions from more than one
   file. `python3 scripts/check_theorems_sigs.py docs/THEOREMS.md` lists the blocks that are
   not verbatim, so the auditor knows which ones to read next to the source: 148 signatures,
   136 verbatim, 12 no verbatim match. The script skips a block whose signature is shorter
   than 25 characters or a name it already compared (`scripts/check_theorems_sigs.py:23`),
   so its total can be lower than the number of `theorem`/`def`/`structure` blocks quoted in
   this document. All twelve are accounted for. Ten are model structures or definitions in
   section 1, where only the implicit size and space binders are dropped (the preamble of
   section 1 names them). Two are the Gaussian blocks of section 6.4, abridged inside the
   conclusion at the elaborated proof terms; the "Hypotheses in words" there spells out what
   each `…` and `_` stands for and lists every binder. **No binder is dropped from the
   source text of any theorem signature in this document.** The instance `[NeZero M]` (at
   least one table; `MultiTableModel.stack` reads `(m.tbl 0).v`, `StackSVD.lean:157`) is
   declared once per section by `variable [NeZero M]` in 13 rank-one files
   (`StackSVD.lean:74`, `StackSVDWeighted.lean:1064`, `SVDStack/Main.lean:95`,
   `SVDStack/Defs.lean:476`, `SVDStack/Weighted.lean:528, 685`, `SVDStack/Gram.lean:429`,
   `SVDStack/Rayleigh.lean:313`, `SVDStack/Simple.lean:50`, `MLE.lean:115`,
   `MLEConverse.lean:120`, `RMT/Het/EdgeSharp.lean:774`, `RMT/Het/Sup.lean:104`,
   `StrictFacades.lean:40`). Lean adds it to every declaration of the section whose
   statement mentions `M`, whether or not that declaration repeats it. For 26 rank-one
   declarations (listed in `notes/archive/F22_nezero.md`) the quoted block also writes
   `[NeZero M]` explicitly, as its first binder, with `omit [NeZero M] in` in front so the
   section instance is not duplicated; that way the signature quoted here, read on its own,
   already shows every instance its elaborated type carries. `#check` on the built tree
   confirms exactly one `[NeZero M]` per declaration. Nothing else is hidden: no other
   section-level instance or variable enters a quoted declaration without appearing in its
   source text. The hypothesis is well-formedness: at least one table.
4. **Hypotheses in words.** One line per binder.
5. **Modeling choices.** Numbered, from the archived review note and from `notes/FLAGGED.md`.
6. **Proof sketch.** Where present, the route of the proof, with the lemmas it passes through
   and their file lines.
7. **Scope.** What the Lean statement says that the paper's does not, or the reverse.
8. **Audit trail.** The audit notes that examined the result, and their verdicts.

Two forms exist for almost every result.

- **Layer 1** takes the random matrix theory input as an explicit hypothesis structure
  (`SingleTableLaw`, `HeteroLaw`, `TableLawR`, `SubspaceLaw`, `SubspaceLawG`, `HeteroLawR`).
  The theorem is an honest implication and holds for any noise law that satisfies the
  structure.
- **Gaussian**, named `*_gaussian`, takes model hypotheses only: independent Gaussian tables
  in the proportional regime. The structure is proved, not assumed.

A third form exists for five results since 2026-09-09. **General noise**, named `*_of_general`,
takes model hypotheses at the paper's own noise class `assum:general_noise`: independent tables
with i.i.d. entries of one fixed law with mean `0`, variance `1` and a finite fourth moment.
It adds two binders that the paper does not have, a Lebesgue density and the upper edge of the
bulk, so it is an implication and not an unconditional theorem. Since 2026-09-10 evening
(Stage 3 of `notes/NONGAUSSIAN_SCOPE.md`) the tree proves the edge binder at four moments
unconditionally (`SpikedModel.lamMax_W0_edge_of_general`, `RMT/General/Edge/Sup.lean`;
section 7.1 names the one combinatorial count the proof needed, proved 2026-09-10 evening),
and a fourth form,
named `*_of_moments`, takes the same binders minus the edge. Only the density binder remains
beyond the paper. Sections 2.1 and 3.1 give both binders, and `notes/NONGAUSSIAN_SCOPE.md`
gives the staged plan they come from.

There is no `axiom` in the development. Convergence is always in probability, never almost
surely, and that is the paper's own mode: `main_paper.tex:234` defines `\pto` as "convergence
in probability of a sequence of random variables", and `:159` and `:175` state the results in
that mode. The paper writes "almost surely" twice (`:1129`, `:1142`), both inside one proof
and not in any statement. The formalization therefore does not weaken the paper's claim.

### Glossary of the project's own terms

An auditor meets these words below. None of them is standard Lean or standard random matrix
theory; they name parts of this development.

| Term | Meaning |
|---|---|
| **Layer 1** | A theorem that takes the random matrix theory input as an explicit hypothesis structure. It is an honest implication and assumes no noise law. |
| **Layer 2** | The Gaussian proof of such a structure, that is a `*_of_gaussian` theorem. Layer 1 plus Layer 2 gives a theorem with model hypotheses only. |
| **facade** | A theorem that restates a proved result in the paper's own shape (inner product in place of projector, one weighting in place of another, a conjunction in place of two lemmas). It adds no mathematics. Facade names end in `_inner`, `_paper`, `_zero`, `_opt`, `_full`, `_frobenius`, `_eig`. |
| **`*_of_gaussian`** | The name of a Layer 2 theorem: it proves a hypothesis structure for independent Gaussian tables in the proportional regime. |
| **`*_gaussian`** | The name of a result that takes model hypotheses only, obtained by feeding a `*_of_gaussian` theorem to a Layer 1 theorem. |
| **`*_of_general`** | A general-noise name: a structure discharged at one fixed i.i.d. law (`thetaEstLaw_of_general`, Stage 0; `singleTableLaw_of_general`, Stage 1), or a Layer 1 theorem fed by that discharge (the four facades of `General/Layer1.lean`). Every Stage 1 name carries the upper edge of the bulk as a hypothesis. |
| **`*_of_moments`** | The same result with the edge hypothesis supplied (Stage 3, 2026-09-10). It takes the paper's noise class and a Lebesgue density, and nothing else, but it rests on one open combinatorial count (section 7.1). Five names: `singleTableLaw_of_moments` (`RMT/General/Edge/Sup.lean`) and the four facades of `General/Layer1.lean`. |
| **Stage 0 to 5** | The stages of `notes/NONGAUSSIAN_SCOPE.md`, the plan that carries the rank-one results from Gaussian noise to the paper's four-moment class. Stage 0 (`ThetaEstLaw`), Stage 1 (`SingleTableLaw`, given the edge) and Stage 3 (the edge itself, at four moments) are done. |
| **`rk`** | The per-table spike count `r_i` of `assum:unaligned`, as a function `Fin M → ℕ`. `alignedRk M r` is the constant function `r`, that is the exactly aligned family of `assum:rank_r`. |
| **item R0 to R6, S, Sym, T** | Stages of the Marchenko-Pastur chain that proves `singleTableLaw_of_gaussian`, in `RMT/`. The roadmap is `notes/archive/rmt_roadmap.md`. |
| **item P, item D** | Two roadmap items outside that chain: P is `lem:noise_projection_concentration` (`ThetaEst.lean`), D is `lem:delocalization` (`SVDStack/Gram.lean`). |
| **H0 to H12** | The same for the heteroscedastic chain `RMT/Het/` that proves `heteroLaw_of_gaussian`. |
| **Track A to E** | The five work streams of Section 7: A (Gaussian discharge of `SubspaceLaw` at `r_i = 1`), B (Section 7 at general `r_i`, Layer 1 only), C (Gaussian discharge of `TableLawR` at general `r_i`), D (the aligned rank-`r` svdstack clauses), E (Gaussian discharge of `HeteroLawR`). A "wave" is one batch of agents inside a track. |
| **D1 to D38** | Numbered decisions Claude took alone, in `notes/FLAGGED.md`. D32 carries a numbered sub-list, cited here as "item *n* of the D32 block". |
| **Q1 to Q17** | Numbered questions for the user, in the same file. |
| **E1 to E12** | Numbered findings about the paper, in `notes/paper_edits.md`. Each is a proposal; none is applied to the paper. |
| **class A, B, C** | The formalization class of a result in `notes/README.md`: A scalar algebra, B finite matrix algebra, C a random matrix theory limit law. |

`lean/StackedSVD/Main.lean` (added 2026-09-07) collects the final form of every result
above into one file, in the order of this document. Each paper label gets a theorem named
after it (`prop_single_table`, `lem_delocalization`, ..., `prop_singleweight_suboptimality`
of 6.9), whose signature is the `_gaussian` theorem quoted in that section, or the
deterministic theorem where the section has one instead. Its proof is one call to that
theorem. Two restatements, `lem_delocalization` (3.2) and `lem_general_rank_delocalization`
(6.1), compose the Layer 1 theorem quoted in that section with the `_of_gaussian` theorems
that discharge its `law` hypothesis, because the tree has no single theorem for their
Gaussian form. `Main.lean` carries no other proof, so a reader can see every result of the
paper in its final form from that one file, with the checker guaranteeing that each
statement there is exactly what the tree proves.

## 1. Definitions

Throughout this section the **implicit** binders that carry the ambient sizes and spaces are
dropped from the code blocks, and only from them: `{Ω : ℕ → Type*}`,
`[∀ N, MeasurableSpace (Ω N)]`, `{μ : ∀ N, Measure (Ω N)}`, `{M : ℕ}`, `{n d : ℕ}`,
`{r rk : ℕ}`. Lean infers each of them from a later explicit argument. Nothing else is
elided; `python3 scripts/check_theorems_sigs.py docs/THEOREMS.md` names the ten blocks below where
this happens, plus the two abridged blocks of section 6.4.

### 1.1 The models

```lean
/-- Rank-one spiked model on one probability space per `N`: `X N = θ u_N v_Nᵀ + d_N^{-1/2} Z N`.
`Z` is the unscaled noise (entries `N(0,1)` under `GaussianNoise`). -/
structure SpikedModel {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    (μ : ∀ N, Measure (Ω N)) (n d : ℕ → ℕ) where
  θ : ℝ
  u : (N : ℕ) → EuclideanSpace ℝ (Fin (n N))
  v : (N : ℕ) → EuclideanSpace ℝ (Fin (d N))
  Z : (N : ℕ) → Ω N → Matrix (Fin (n N)) (Fin (d N)) ℝ
  hθ : 0 ≤ θ
  hn : ∀ N, 0 < n N
  hd : ∀ N, 0 < d N
  hu : ∀ N, ‖u N‖ = 1
  hv : ∀ N, ‖v N‖ = 1
  hZ : ∀ N, Measurable (Z N)
```

This is `eq:model_rank1` with one asymptotic index `N`. The paper's `n_i, d → ∞` with
`n_i/d → c_i` becomes the predicate `SpikedModel.Regime c`, a conjunction of **three** limits:
`Tendsto n atTop atTop`, `Tendsto d atTop atTop` and
`Tendsto (fun N => (n N : ℝ) / d N) atTop (nhds c)` (`Defs.lean`). `SpikedModelR.Regime`
(`RankR/General.lean`) is the same conjunction for a rank-`rk` table.

```lean
/-- `M` tables with a shared `v`, all on `(Ω N, μ N)`. Cross-table independence is a
separate predicate `JointGaussianNoise`. -/
structure MultiTableModel {Ω : ℕ → Type*} [∀ N, MeasurableSpace (Ω N)]
    (μ : ∀ N, Measure (Ω N)) (M : ℕ) (n : Fin M → ℕ → ℕ) (d : ℕ → ℕ) where
  tbl : (i : Fin M) → SpikedModel μ (n i) d
  hv : ∀ i j N, (tbl i).v N = (tbl j).v N
```

`MultiTableModel.JointGaussianNoise` says the joint law of the `M` noise matrices at level
`N` is the product of `M` standard Gaussian matrix laws. This is `assum:general_noise` in
the Gaussian case, plus the cross-table independence that the paper states in words and that
a per-table predicate cannot express.

```lean
/-- One table of `assum:unaligned` with `rk` spikes: `X = U Θ Vᵀ + d^{-1/2} Z` with `U` and `V`
of orthonormal columns and `Θ = diag(θ_1, …, θ_rk)`. The rank-`rk` twin of `SpikedModel`. -/
structure SpikedModelR (μ : ∀ N, Measure (Ω N)) (n d : ℕ → ℕ) (rk : ℕ) where
  θ : Fin rk → ℝ
  U : (N : ℕ) → Matrix (Fin (n N)) (Fin rk) ℝ
  V : (N : ℕ) → Matrix (Fin (d N)) (Fin rk) ℝ
  Z : (N : ℕ) → Ω N → Matrix (Fin (n N)) (Fin (d N)) ℝ
  hθnn : ∀ k, 0 ≤ θ k
  hθanti : StrictAnti θ
  hn : ∀ N, 0 < n N
  hd : ∀ N, 0 < d N
  hU : ∀ N, (U N)ᵀ * U N = 1
  hV : ∀ N, (V N)ᵀ * V N = 1
  hZ : ∀ N, Measurable (Z N)
```

`hθnn : ∀ k, 0 ≤ θ k` and `hθanti : StrictAnti θ` fix the spikes of a table: nonnegative and
strictly decreasing. So every spike is positive except possibly the last one, which may be
`0`. A table whose last spike is `0` is the same table with one spike fewer, so this class
is the paper's class of `assum:unaligned` ("diagonal with positive entries") up to that
relabeling; the nonnegative form makes the rank-one table `SpikedModel` (whose `hθ : 0 ≤ θ`
allows `θ = 0`) a special case through `SpikedModel.toRankR`. `hθanti` is the paper's
distinctness hypothesis (`main_paper.tex:768`) together with a fixed column order, which is
free: a simultaneous permutation of the columns of `U_i`, `Θ_i` and `R_i` leaves the model
unchanged. Both are model fields, not hypotheses of any theorem, and no Layer 1 proof reads
either; see section 6.7 and `notes/paper_edits.md` E7.

```lean
/-- `assum:unaligned` with no restriction on the `r_i`. -/
structure UnalignedModelR (μ : ∀ N, Measure (Ω N)) (M : ℕ) (n : Fin M → ℕ → ℕ) (d : ℕ → ℕ)
    (r : ℕ) (rk : Fin M → ℕ) where
  tbl : (i : Fin M) → SpikedModelR μ (n i) d (rk i)
  V : (N : ℕ) → Matrix (Fin (d N)) (Fin r) ℝ
  R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ
  hV : ∀ N, (V N)ᵀ * V N = 1
  hR : ∀ i, (R i)ᵀ * R i = 1
  hv : ∀ i N, (tbl i).V N = V N * R i
```

`UnalignedModel μ M n d r` (`RankR/Defs.lean`) is the `r_i = 1` model. It is a **separate**
structure, built on `SpikedModel` tables, with `R : Fin M → EuclideanSpace ℝ (Fin r)` a unit
vector per table in place of a matrix; the fields `hθnn` and `hθanti` do not occur in it. It
is not `UnalignedModelR` at `rk = 1`. The tree carries one direction of a bridge:
`UnalignedModelR.toUnaligned` reads a model at `rk = fun _ => 1` as an `UnalignedModel`, and
`perfRG_one_eq_perfR` identifies the two performance measures (`RankR/GeneralFrob.lean:632`,
`:701`); the table-level map `SpikedModel.toRankR` (`RankR/General.lean`) goes the other way
for one table, but no model-level map `UnalignedModel → UnalignedModelR` exists. Both twins
are proved separately, and both are listed in every section below. The paper's `Rank(∑_i R_i R_iᵀ) = r`
is **not** a field: the unweighted proposition does not use it, and the weighted theorem
takes the stronger `rank B_R = r` as an explicit argument. The weaker pair "`∃ k, 0 < β k`
and `Rank(∑_i R_i R_iᵀ) = r`" does not suffice (counterexample in `notes/FLAGGED.md`, D16).

`alignedRk M r` is the constant map `fun _ => r`. The exactly aligned family of
`assum:rank_r` is `UnalignedModelR μ M n d r (alignedRk M r)` with `R_i = 1`.

### 1.2 Convergence and spectral objects

```lean
/-- Convergence in probability to a constant along `N → ∞`. -/
def TendstoInProb (μ : ∀ N, Measure (Ω N)) (f : ∀ N, Ω N → ℝ) (a : ℝ) : Prop :=
  ∀ ε > 0, Tendsto (fun N => μ N {ω | ε ≤ |f N ω - a|}) atTop (𝓝 0)

/-- Squared overlap of the top right singular subspace of `X` with `w`. -/
noncomputable def overlap (X : Matrix (Fin n) (Fin d) ℝ) (w : EuclideanSpace ℝ (Fin d)) : ℝ :=
  ‖topProj (Xᵀ * X) (isHermitian_transpose_mul_self X) w‖ ^ 2
```

`overlap` is the projector form. It avoids the sign ambiguity of an eigenvector. On the event
where the top eigenvalue is simple it equals the paper's `⟨v̂, w⟩²`
(`overlap_eq_inner_sq`). The headline estimator theorems (`prop:stacksvd_general`,
`thm:svd_stack_general`, `thm:stacksvd_weighted`, `thm:svdstack_weighted`, the rank-`r`
results) carry an inner-product twin named `_inner` that consumes that fact; the
specializations that are stated in projector form only (`prop:single_table`, `thm:simple_thm1`,
`cor.2`, the inadmissibility and remark instances, the theta estimator) are one application of
`overlap_eq_inner_sq` away from the paper's `⟨v̂, w⟩²`, and the tree does not repeat that
step for them (external rank-one audit finding 8). `overlapIdx X k w` is the same quantity at the `k`-th eigenvalue,
used by the rank-`r` results. `SimpleIdx A hA k` says eigenvalue `k` is simple;
`SimpleSpec A hA r` says the top `r` are.

### 1.3 The scalars

```lean
/-- `β²` of `prop:single_table`. -/
noncomputable def betaSq (θ c : ℝ) : ℝ :=
  if θ ^ 4 > c then (θ ^ 4 - c) / (θ ^ 4 + θ ^ 2) else 0

noncomputable def beta (θ c : ℝ) : ℝ := Real.sqrt (betaSq θ c)

/-- Bulk edge of `W = Eᵀ E` under the `1/d` variance scaling. -/
noncomputable def bulkEdge (c : ℝ) : ℝ := (1 + Real.sqrt c) ^ 2

/-- Limit of the top eigenvalue of `Xᵀ X`. -/
noncomputable def rhoSq (θ c : ℝ) : ℝ :=
  if θ ^ 4 > c then θ ^ 2 + 1 + c + c / θ ^ 2 else bulkEdge c

/-- `A_β = β βᵀ + diag(1 - β_1², ..., 1 - β_M²)` (`eq:A_beta_main_text`). -/
noncomputable def Abeta (β : Fin M → ℝ) : Matrix (Fin M) (Fin M) ℝ :=
  Matrix.vecMulVec β β + Matrix.diagonal fun i => 1 - β i ^ 2

/-- Closed-form limit of unweighted stacksvd (`prop:stacksvd_general`). -/
noncomputable def stackSVDLimit (θ c : Fin M → ℝ) : ℝ :=
  if (∑ i, θ i ^ 2) ^ 2 > ∑ i, c i then
    ((∑ i, θ i ^ 2) ^ 2 - ∑ i, c i) / ((∑ i, θ i ^ 2) * ((∑ i, θ i ^ 2) + 1))
  else 0

/-- Limit value of `thm:svd_stack_general`: `(βᵀ v_max(A_β))² / λ_max(A_β)`. -/
noncomputable def svdstackLimit (β : Fin M → ℝ) : ℝ :=
  (β ⬝ᵥ WithLp.ofLp (vMax (Abeta β) (isHermitian_Abeta β))) ^ 2
    / lamMax (Abeta β) (isHermitian_Abeta β)

/-- `S = ∑_i β_i² / (1 - β_i²)` (`thm:svdstack_weighted`). -/
noncomputable def Sval (β : Fin M → ℝ) : ℝ := ∑ i, β i ^ 2 / (1 - β i ^ 2)

/-- The optimal svdstack weights `w_i⋆ = 1/√(1 - β_i²)`. -/
noncomputable def optW (β : Fin M → ℝ) : Fin M → ℝ := fun i => 1 / Real.sqrt (1 - β i ^ 2)

/-- The limit of optimally weighted svdstack, `S/(S+1)`. -/
noncomputable def svdstackLimitOpt (β : Fin M → ℝ) : ℝ := Sval β / (Sval β + 1)
```

The threshold convention differs from the paper on purpose. The paper writes
`θ ≥ c^{1/4}`; Lean writes `c < θ⁴`. Both branches give `0` at equality, so the two agree;
the equivalence of the two conditions uses the model field `SpikedModel.hθ : 0 ≤ θ`.

Weighted stackSVD scalars live in the `Scalars` namespace: `optWstack θ c` is
`θ_i/√(θ_i²+c_i)`; `Lw θ c w` is the paper's `L(w)`, total, with `0` below `eq:assumption4`;
`stackSVDLimitW θ c` is the unique root of `gW θ c = 1` in `(0,1)` (`Scalars.lean:408`), the
paper's `γ_opt`; `Scalars.L_optW_eq` (`StackSVDWeighted.lean:1003`) identifies it with `Lw` at
the optimal weights;
`Assumption4 θ c w` is `eq:assumption4`; `binaryStackSVDLimit S θ c` is the `cor.2` value on
a subset `S`.

Rank-`r` scalars, all in `Scalars`:

```lean
/-- `w_ij = θ_ij/√(θ_ij² + c_i)` (`eq:stacksvd_app_wij`). -/
noncomputable def wStackR (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (j : Fin r) : Fin M → ℝ :=
  optWstack (fun i => θ i j) c

/-- `θ̃_jk²` of `eq:stacksvd_apptildTheta`. -/
noncomputable def thetaTildeSq (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (j k : Fin r) : ℝ :=
  ∑ i, θ i j ^ 2 * θ i k ^ 2 / (θ i j ^ 2 + c i)

/-- `γ_j` (`eq:stacksvd_gammak`), total: the rank-1 limit `stackSVDLimitW` on column `j`. -/
noncomputable def gammaR (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (j : Fin r) : ℝ :=
  stackSVDLimitW (fun i => θ i j) c

/-- `ℓ_j` of `alg:rank_r_stacksvd`, 0-based. -/
noncomputable def ellR (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (j : Fin r) : ℕ :=
  (Finset.univ.filter fun k => thetaTildeSq θ c j j < thetaTildeSq θ c j k).card

/-- The number of components supercritical at the weights `w`. -/
noncomputable def numSup (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (w : Fin M → ℝ) : ℕ :=
  (Finset.univ.filter fun l : Fin r => Assumption4 (fun i => θ i l) c w).card

/-- The sorted index of the outlier of component `k` at the weights `w`. -/
noncomputable def ellSup (θ : Fin M → Fin r → ℝ) (c : Fin M → ℝ) (w : Fin M → ℝ) (k : Fin r) :
    ℕ :=
  (Finset.univ.filter fun l : Fin r =>
    Assumption4 (fun i => θ i l) c w ∧
      MPhet.rhoHet (fun i => θ i k) c w < MPhet.rhoHet (fun i => θ i l) c w).card

/-- The strength table `θ_ij` of the exactly aligned rank-`r` model. -/
noncomputable def UnalignedModelR.thetaAligned (m : UnalignedModelR μ M n d r (alignedRk M r)) :
    Fin M → Fin r → ℝ := fun i j => (m.tbl i).θ j
```

`ellSup` and `ellR` differ in general. `Scalars.ellSup_eq_ellR` proves they agree inside the
model class of `assum:rank_r`. The difference is paper finding E7; see section 7.

### 1.4 The estimators

| Object | Lean name | File | Meaning |
|---|---|---|---|
| `v̂_stacksvd` (unit weights) | `MultiTableModel.stack`, then `overlap (m.stack.X N ω) (m.stack.v N)` | `StackSVD.lean` | performance of the unweighted stack |
| `v̂_stacksvd(w)` | `MultiTableModel.stackPerfW w` | `StackSVDWeighted.lean` | performance of the weighted stack |
| `v̂_svdstack` | `MultiTableModel.svdstackPerf` | `SVDStack/Defs.lean` | top right singular vector of the stacked `v̂_iᵀ` |
| `v̂_svdstack(w)` | `MultiTableModel.svdstackPerfW w` | `SVDStack/Defs.lean` | weighted svdstack performance |
| `θ̂_2` | `MultiTableModel.thetaHat i j ci cj` | `ThetaEst.lean` | `eq:theta_estimation` |
| `‖V̂_svdstackᵀ V‖_F²`, rank `r` | `UnalignedModelR.perfRG`, `perfRGW W` | `RankR/General.lean` | unweighted and weighted rank-`r` svdstack |
| `‖V̂_stacksvdᵀ V‖_F²`, rank `r` | `UnalignedModelR.perfStackRG` | `RankR/SubspaceG.lean` | rank-`r` stacksvd subspace performance (projector-sum form); the paper-literal `frobSq ((topEigMat ...)ᵀ * V)` is the conclusion of `prop_stacksvd_subspace_general_frobenius_eig_gaussian` (`RankR/SubspaceGaussian.lean`) |
| `v̂_{j,stacksvd}` | `UnalignedModelR.vhatStackR c j` | `RankR/StackGamma.lean` | component `j` of rank-`r` weighted stacksvd |
| `‖Vᵀ V̂_stacksvd‖_F²` | `UnalignedModelR.frobSqStackR c` | `RankR/StackGamma.lean` | `∑_j ∑_k ⟪v̂_j, v_k⟫²` |
| `V̂_svdstack(W)` frame | `UnalignedModelR.vhatSvdstackGW W N ω frame val` | `RankR/GeneralFrob.lean` | the matrix whose Frobenius norm is `perfRGW` |

## 2. Hypothesis structures and their Gaussian discharge

The paper cites three random matrix theory results as black boxes: `liu2023asymptotic`
Theorems 1 and 2, and `10.3150/19-BEJ1129` Theorem 2.3. No prover reaches those today. They
enter this development as `structure`s, one field per fact a proof reads, and each structure
is then proved for Gaussian noise by a named theorem. Nothing is an `axiom`.

### 2.1 `SpikedModel.SingleTableLaw` (`RMT.lean`)

```lean
structure SingleTableLaw (m : SpikedModel μ n d) (c : ℝ) : Prop where
  align : TendstoInProb μ (fun N ω => overlap (m.X N ω) (m.v N)) (betaSq m.θ c)
  delocUniform : ∀ ε > 0,
    Tendsto (fun N => ⨆ w ∈ m.orthUnit N, μ N {ω | ε ≤ overlap (m.X N ω) w}) atTop (𝓝 0)
  lamMax : TendstoInProb μ (fun N ω => gramLamMax (m.X N ω)) (rhoSq m.θ c)
  topSimple : ∀ N, ∀ᵐ ω ∂(μ N),
    TopSimple ((m.X N ω)ᵀ * m.X N ω) (isHermitian_transpose_mul_self (m.X N ω))
```

Fields in words. `align` is part 1 of `prop:single_table`. `delocUniform` is part 2,
strengthened to be uniform over unit `w ⊥ v`; the paper's sequential form is the theorem
`SingleTableLaw.deloc_seq`. `lamMax` is `σ_1²(X) → θ²+1+c+c/θ²`, cited by the paper in
`thm:theta_est`. `topSimple` is finite-`N` almost-sure simplicity of the top eigenvalue,
which the paper does not state; it replaces the paper's implicit "the top singular vector is
well defined".

Discharged by:

```lean
theorem singleTableLaw_of_gaussian [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) (hG : m.GaussianNoise) :
    m.SingleTableLaw c
```

in `RMT/Full.lean`, both regimes (`c < θ⁴` and `θ⁴ ≤ c`), through the chain
`RMT/{MP, MP7, R0, R4, R4C, Simplicity, Symmetry, R3, R3minus, T, R1, ResolvDeriv, SteinStep,
R2, R5, R6, Sup, TailShift}`. The literature mirror is `liu2023asymptotic` Theorem 1.

And, since 2026-09-09, at a general law with a density, given the edge:

```lean
theorem singleTableLaw_of_general [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.GeneralNoise ν)
    (hedge : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ bulkEdge c + ε}) atTop (𝓝 1)) :
    m.SingleTableLaw c
```

in `RMT/General/Sup.lean` (Stage 1 of `notes/NONGAUSSIAN_SCOPE.md`), through the chain
`RMT/General/{Defs, QuadForm, Stability, Companion, MPtilde, Trace, R3minus, Simplicity, Iso,
IsoMixed, FormsBridge, ProbC, Forms, FormsLimits, FormsSup, FormsGeneral, Deloc, DelocLimits,
DelocAlign, DelocUniform, DelocUniformSup, Sup}` and the probability files
`Prob/{NoiseLaw, NoiseMoments, LinFormMoments}` (Stage 0) and `Prob/{Chebyshev, Tensorization}`
(Stage 1). `hν : NoiseLaw ν` is the paper's
four-moment class at one fixed law (section 2.6). `hG : m.GeneralNoise ν` says that the
entries of `Z N` are i.i.d. of law `ν` at every `N`. Two binders go beyond the paper.
`hedge` is the upper edge of the bulk: for every `ε > 0` the probability that the top
eigenvalue of the downdated Gram matrix `W₀` is at most `bulkEdge c + ε` tends to `1`. Since
2026-09-10 evening the tree proves it at four moments unconditionally:
`lamMax_W0_edge_of_general` in
`RMT/General/Edge/Sup.lean` (Stage 3, route C of `notes/stage3_edge.md`: the moment method
with a truncation at `d^(1/2 - 1/16)`, a Markov bound at the order `⌈C log d⌉`, and the
comparison of the truncated trace with the Gaussian one) through the chain
`RMT/General/Edge/{Defs, Arith, Trunc, Sparse, Trace, Compare, Gaussian, Count, Code, Dyck,
CountBound, Excess, Markov, Sup}`, including the one combinatorial count that section 7.1
names, proved 2026-09-10 evening. The theorem with `hedge` stays as the Stage 1 statement; the closed form is

```lean
theorem singleTableLaw_of_moments [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.GeneralNoise ν) :
    m.SingleTableLaw c
```

in the same file. The Gaussian chain gets the edge from Sudakov-Fernique (`RMT/R3.lean`), a
tool no general law admits; the general chain gets it from the trace of `(ZᵀZ)^k`.
`hac : ν ≪ volume` says that `ν` has a Lebesgue density. It pays for `topSimple`, which is
false for Rademacher entries, and for one step of the subcritical `align`. Section 3.1 gives
the detail and the route.

### 2.2 `MultiTableModel.HeteroLaw` (`StackSVDWeighted.lean`)

```lean
structure HeteroLaw [NeZero M] (m : MultiTableModel μ M n d) (w c : Fin M → ℝ) : Prop where
  align : TendstoInProb μ (fun N ω => m.stackPerfW w N ω)
    (Scalars.Lw (fun i => (m.tbl i).θ) c w)
  topSimple : ∀ N, ∀ᵐ ω ∂(μ N),
    TopSimple (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω)
```

`align` is the specialization of `eq:weighted_norm` at `b = v`, with both branches of the
paper inside `Lw`. The paper's `b`-general half is not a field: no result below consumes it.
The literature mirror is `10.3150/19-BEJ1129` Theorem 2.3.

Discharged by:

```lean
theorem heteroLaw_of_gaussian [NeZero M]
    (m : MultiTableModel μ M n d) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    m.HeteroLaw w c
```

in `RMT/Het/Sup.lean`. It takes no condition beyond `hc`, `hreg`, `hG` and a nonzero weight
vector `hw : ∃ i, w i ≠ 0`; no margin or separation condition. The exact heteroscedastic upper edge, which two
expert reviews put out of reach through Bai-Silverstein, is proved separately by a sharp
Sudakov-Fernique bound:

```lean
theorem heteroEdge_of_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hw : ∃ i, w i ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    m.HeteroEdge w c (MPhet.bHet c w)
```

in `RMT/Het/EdgeSharp.lean`. `HeteroEdge` stays in the tree as the general interface for
another noise law. `heteroLaw_of_gaussian_margin` is the weaker interim theorem, kept as a
record; it needs `MPhet.bSF c w < MPhet.rhoHet θ c w`.

### 2.3 `SpikedModelR.TableLawR` (`RankR/General.lean`)

```lean
structure TableLawR (m : SpikedModelR μ n d rk) (c : ℝ) : Prop where
  align : ∀ k : Fin rk,
    TendstoInProb μ (fun N ω => overlapIdx (m.X N ω) (k : ℕ) (m.col N k)) (betaSq (m.θ k) c)
  cross : ∀ k l : Fin rk, k ≠ l →
    TendstoInProb μ (fun N ω => overlapIdx (m.X N ω) (k : ℕ) (m.col N l)) 0
  delocUniform : ∀ k : Fin rk, ∀ ε > 0,
    Tendsto (fun N => ⨆ w ∈ m.orthUnitR N, μ N {ω | ε ≤ overlapIdx (m.X N ω) (k : ℕ) w})
      atTop (𝓝 0)
  simple : ∀ N, ∀ᵐ ω ∂(μ N),
    SimpleSpec ((m.X N ω)ᵀ * m.X N ω) (isHermitian_transpose_mul_self (m.X N ω)) rk
```

`prop:single_table` at rank `rk`, in projector form. `cross` is where the distinctness of the
`θ_k` is used; it is false at a tie. There is no `lamMax` twin, because no statement of
Section 7 reads it.

Discharged by:

```lean
theorem SpikedModelR.tableLawR_of_gaussian_rk [∀ N, IsProbabilityMeasure (μ N)]
    (m : SpikedModelR μ n d rk) {c : ℝ} (hc : 0 < c) (hreg : m.Regime c)
    (hG : m.GaussianNoise) : m.TableLawR c

theorem UnalignedModelR.tableLawR_of_gaussian_rk [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) {c : Fin M → ℝ} (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    ∀ i, (m.tbl i).TableLawR (c i)
```

both in `RankR/RMT/TableLawGaussian.lean`. No supercriticality is assumed; the proof case
splits per index on `c < θ_k⁴`.

### 2.4 `UnalignedModel.SubspaceLaw` and `UnalignedModelR.SubspaceLawG`

```lean
structure SubspaceLaw (m : UnalignedModel μ M n d r) (c : ℝ) : Prop where
  align : ∀ j : Fin r, TendstoInProb μ
    (fun N ω => ‖specProjTop (m.stackGram N ω) (m.isHermitian_stackGram N ω) r
      (m.spikeVec j N)‖ ^ 2) (betaSq (Real.sqrt (m.coreEig j)) c)
```

The rank-`r` spiked law of the unit-weight stack, at aspect ratio `‖c‖₁ = ∑_i c_i`. It is
`prop:single_table` read once per spike of the core matrix `C = ∑_i R_i Θ_i² R_iᵀ`. The
structure carries no eigengap of `C` and no distinct-spike condition (decisions D17 and D18).
`SubspaceLawG` (`RankR/SubspaceG.lean`) is the same at general `r_i`.

Discharged by `UnalignedModel.subspaceLaw_of_gaussian` and
`UnalignedModelR.subspaceLawG_of_gaussian`, both in `RankR/SubspaceGaussian.lean`:

```lean
theorem subspaceLaw_of_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ M n d r) (hG : m.JointGaussianNoise)
    {cc : Fin M → ℝ} (hreg : ∀ i, (m.tbl i).Regime (cc i)) (hc : 0 < ∑ i, cc i) :
    m.SubspaceLaw (∑ i, cc i)
```

There is no side condition on `∑_i n_i`: the tail shift of `RankR/RMT/ShiftR.lean` removes it.

### 2.5 `UnalignedModelR.HeteroLawR` (`RankR/StackGamma.lean` since F28, 2026-09-08; was `RankR/StackMain.lean`, which still holds its consumer theorems)

```lean
structure HeteroLawR (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ) :
    Prop where
  align : ∀ j : Fin r, TendstoInProb μ (fun N ω => m.stackOverlapJ c j N ω)
    (Scalars.gammaR m.thetaAligned c j)
  crossProj : ∀ j k : Fin r, k ≠ j →
    TendstoInProb μ (fun N ω => overlapIdx (m.stackXJ c j N ω)
      (Scalars.ellR m.thetaAligned c j) (m.colVecG N k)) 0
  simpleIdxJ : ∀ (j : Fin r) (N : ℕ), ∀ᵐ ω ∂(μ N),
    SimpleIdx (m.stackGramJ c j N ω) (m.isHermitian_stackGramJ c j N ω)
      (Scalars.ellR m.thetaAligned c j)
```

The rank-`r` twin of `HeteroLaw`, one field per thing `thm:rank_r_stacksvd` reads. Every
field is a statement about the `r` weighted stacks `X_stack^{(j)}` of
`eq:stacksvd_appXstack`. `align` is the first display of the corollary in projector form.
`crossProj` is the paper's sentence "the columns of `V̂` will be asymptotically orthogonal"
(`main_paper.tex:2318`), which the second display needs. `simpleIdxJ` turns the projector
overlap into the paper's `⟪v̂_j, v_j⟫²`.

Discharged by:

```lean
theorem heteroLawR_of_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (hnz : ∀ j : Fin r, ∃ i, 0 < (m.tbl i).θ j) :
    m.HeteroLawR c
```

in `RankR/Het/Sup.lean`, the last stage of Track E. The chain is `RankR/Het/{Scalars, Split,
Duality, Simplicity, Edge, Forms, Deloc, BulkDet, Outliers, Align, Bulk, Sup}`. It uses two
further internal structures: `HeteroEdgeR` (`RankR/Het/Edge.lean`), the rank-`r` edge, proved
at `b = MPhet.bHet c w` by `heteroEdgeR_of_gaussian`; and `ResolventLimitsHetR`
(`RankR/Het/Forms.lean`), the rank-`r` (H1) and (H2) resolvent limits, proved by
`resolventLimitsHetR_of_gaussian`.

### 2.6 `MultiTableModel.ThetaEstLaw` (`ThetaEst.lean`)

```lean
structure ThetaEstLaw (m : MultiTableModel μ M n d) (i j : Fin M) (cj : ℝ) : Prop where
  noiseProj : TendstoInProb μ (fun N ω => ThetaEst.sqNormMulVec ((m.tbl j).E N ω)
    (ThetaEst.topDir ((m.tbl i).X N ω))) cj
  crossTerm : TendstoInProb μ (fun N ω => WithLp.ofLp ((m.tbl j).u N) ⬝ᵥ
    ((m.tbl j).E N ω *ᵥ WithLp.ofLp (ThetaEst.topDir ((m.tbl i).X N ω)))) 0
```

Fields in words (F31, 2026-09-08). Both are limits of the noise `E_j` of table `j` against
the top right singular direction `v̂_i` of table `i`, the two concentration steps that the
paper's proof of `thm:theta_est` takes from `lem:noise_projection_concentration` and the
independence of the tables. `noiseProj` is `‖E_j v̂_i‖² → c_j`, item P of the roadmap with a
random direction in place of a fixed one. `crossTerm` is `u_jᵀ E_j v̂_i → 0`. They are the
only probabilistic input of `thm_theta_est` beyond the `SingleTableLaw` of table `i`.

Discharged by:

```lean
theorem thetaEstLaw_of_gaussian (m : MultiTableModel μ M n d) {i j : Fin M} (hij : i ≠ j)
    (hG : m.JointGaussianNoise) {cj : ℝ} (hregj : (m.tbl j).Regime cj) :
    m.ThetaEstLaw i j cj
```

in `ThetaEst.lean`, from `noise_projection_topDir_tendsto` and `cross_term_tendsto` of the
same file: Chebyshev bounds after the Fubini step `measure_randomDir_le`, where the direction
`v̂_i` is a measurable function of the noise of table `i`, independent of table `j`.

Also discharged for a general i.i.d. noise law (Stage 0 of `notes/NONGAUSSIAN_SCOPE.md`,
2026-09-08). The law is one probability measure `ν` on `ℝ`, the same for every `N`, with mean
`0`, second moment `1` and an integrable fourth power:

```lean
structure NoiseLaw (ν : Measure ℝ) : Prop where
  prob : IsProbabilityMeasure ν
  mean : ∫ x, x ∂ν = 0
  var : ∫ x, x ^ 2 ∂ν = 1
  mom4 : Integrable (fun x => x ^ 4) ν

noncomputable def noiseMatrix (ν : Measure ℝ) (n d : ℕ) : Measure (Matrix (Fin n) (Fin d) ℝ) :=
  Measure.pi fun _ : Fin n => Measure.pi fun _ : Fin d => ν

def MultiTableModel.JointGeneralNoise (m : MultiTableModel μ M n d) (ν : Measure ℝ) : Prop :=
  ∀ N, HasLaw (fun ω i => (m.tbl i).Z N ω)
    (Measure.pi fun i : Fin M => noiseMatrix ν (n i N) (d N)) (μ N)
```

(`Prob/NoiseLaw.lean`; `gaussianMatrix n d = noiseMatrix (gaussianReal 0 1) n d` by `rfl`, and
`m.JointGaussianNoise ↔ m.JointGeneralNoise (gaussianReal 0 1)` by `Iff.rfl`, so the Gaussian
model is the instance `ν = gaussianReal 0 1`, `noiseLaw_gaussian`.) The discharge is
`thetaEstLaw_of_general` (section 5.4), with the same two Chebyshev steps: the Gaussian
exact variance `R2.varSq * n_j / d²` of `‖E_j v̂_i‖²` becomes the bound `(ν₄ + 2) n_j / d²`
from the fourth moment of a linear form (`Prob/LinFormMoments.lean`), and the cross term
keeps the exact second moment `1` (`Prob/NoiseMoments.lean`).

### 2.7 Summary

| Structure | File | Gaussian discharge | Unconditional? |
|---|---|---|---|
| `SpikedModel.SingleTableLaw` | `RMT.lean` | `singleTableLaw_of_gaussian` (`RMT/Full.lean`); also `singleTableLaw_of_general` for any i.i.d. law with mean 0, variance 1, a finite fourth moment and a Lebesgue density (`RMT/General/Sup.lean`, 2026-09-09) | yes for Gaussian; yes at that general law since Stage 3 (`singleTableLaw_of_moments`, `RMT/General/Edge/Sup.lean`, 2026-09-10), which supplies the upper edge of the bulk; the conditional form `singleTableLaw_of_general` keeps the edge as a binder |
| `MultiTableModel.HeteroLaw` | `StackSVDWeighted.lean` | `heteroLaw_of_gaussian` (`RMT/Het/Sup.lean`) | yes, given `hw : ∃ i, w i ≠ 0` |
| `MultiTableModel.HeteroEdge` | `RMT/Het/R4het.lean` | `heteroEdge_of_gaussian` (`RMT/Het/EdgeSharp.lean`) | yes |
| `SpikedModelR.TableLawR` | `RankR/General.lean` | `tableLawR_of_gaussian_rk` (`RankR/RMT/TableLawGaussian.lean`) | yes |
| `UnalignedModel.SubspaceLaw` | `RankR/Subspace.lean` | `subspaceLaw_of_gaussian` (`RankR/SubspaceGaussian.lean`) | yes |
| `UnalignedModelR.SubspaceLawG` | `RankR/SubspaceG.lean` | `subspaceLawG_of_gaussian` (`RankR/SubspaceGaussian.lean`) | yes |
| `UnalignedModelR.HeteroLawR` | `RankR/StackGamma.lean` | `heteroLawR_of_gaussian` (`RankR/Het/Sup.lean`) | yes |
| `UnalignedModelR.HeteroEdgeR` | `RankR/Het/Edge.lean` | `heteroEdgeR_of_gaussian` (`RankR/Het/Edge.lean`) | yes |
| `UnalignedModelR.ResolventLimitsHetR` | `RankR/Het/Forms.lean` | `resolventLimitsHetR_of_gaussian` (`RankR/Het/Forms.lean`) | yes |
| `UnalignedModelR.SingleWeightLaw` | `RankR/SingleWeight/Main.lean` | `singleWeightLaw_of_gaussian` (`RankR/SingleWeight/Het/Sup.lean`) | yes, under the paper's separation assumption `EigSep` (a condition on the parameters, section 6.8) and `hrn` |
| `MultiTableModel.ThetaEstLaw` | `ThetaEst.lean` | `thetaEstLaw_of_gaussian` (`ThetaEst.lean`); also `thetaEstLaw_of_general` for any i.i.d. law with mean 0, variance 1 and a finite fourth moment (`NoiseLaw`, 2026-09-08) | yes |

Every limit-law hypothesis structure in the tree has a Gaussian discharge (`EigSep` of section
6.8 is a condition on the parameters, not a law, and stays a binder). No result of the paper that
this development covers rests on an undischarged black box.

Two structures also have a general-noise discharge, at the paper's own class
`assum:general_noise`: `ThetaEstLaw` unconditionally (`thetaEstLaw_of_general`, Stage 0,
2026-09-08) and `SingleTableLaw` (`singleTableLaw_of_general`, Stage 1, 2026-09-09, which
asks for a Lebesgue density and for the upper edge of the bulk). The edge was a hypothesis
and never a black box: it is a statement about the model, and every theorem that carried it
was an honest implication. Stage 3 (2026-09-10) proves it at four moments, so
`singleTableLaw_of_moments` (`RMT/General/Edge/Sup.lean`) states the same conclusion with the
density as its only binder beyond the paper. One combinatorial input of that proof is stated
and used but not yet proved; section 7.1 names it, and it is the one row of
`docs/SORRIES.md`. Section 3.1 gives the route, and `notes/NONGAUSSIAN_SCOPE.md` gives the
stages.

## 3. Rank-one core results

### 3.1 `prop:single_table`

**Paper.**

```latex
\begin{prop}\label{prop:single_table}
    Under \Cref{assum:general_noise,assum:main}, we have
    \begin{equation*}
        |\langle \hat{v}_i, v \rangle|^2 \ip \beta_i^2 :=  \begin{cases}
        \frac{\theta_i^4 - c_i}{\theta_i^4 + \theta_i^2} & \text{ if } \theta_i \geq c_i^{1/4}, \\
        0 & \text{ otherwise.}
        \end{cases}
    \end{equation*}
\noindent Furthermore, for any deterministic sequence of unit vectors $w^{(d)}$ orthogonal to $v$:
    \begin{equation*}
        |\langle \hat{v}_i, w^{(d)} \rangle|^2 \ip 0.
    \end{equation*}
\end{prop}
```

**Lean.** The result is the hypothesis structure `SingleTableLaw` of section 2.1, and its
Gaussian proof:

```lean
/-- Layer 2 target: the single-table law holds under Gaussian noise in the proportional
regime. -/
theorem singleTableLaw_of_gaussian [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) (hG : m.GaussianNoise) :
    m.SingleTableLaw c
```

**Lean, general noise** (`RMT/General/Sup.lean`, Stage 1 of `notes/NONGAUSSIAN_SCOPE.md`,
2026-09-09; the definitions `NoiseLaw` and `noiseMatrix` are quoted in section 2.6):

```lean
theorem singleTableLaw_of_general [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.GeneralNoise ν)
    (hedge : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ bulkEdge c + ε}) atTop (𝓝 1)) :
    m.SingleTableLaw c
```

**Lean, the edge at four moments** (`RMT/General/Edge/Sup.lean`, Stage 3 of
`notes/NONGAUSSIAN_SCOPE.md`, 2026-09-10). `opNorm_sq_edge_of_general` is the upper half of
the Bai-Yin theorem for the unscaled noise matrix, `lamMax_W0_edge_of_general` is the binder
`hedge` above word for word, and the form without `hedge`, `singleTableLaw_of_moments`, has
the same binders minus the edge; it is `singleTableLaw_of_general` applied to
`lamMax_W0_edge_of_general`:

```lean
theorem opNorm_sq_edge_of_general [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hG : m.GeneralNoise ν) :
    ∀ ε > 0, Tendsto (fun N => μ N
      {ω | ‖m.Z N ω‖ ^ 2 / (d N : ℝ) ≤ bulkEdge c + ε}) atTop (𝓝 1)

theorem lamMax_W0_edge_of_general [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hG : m.GeneralNoise ν) :
    ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (m.W0 N ω) (m.isHermitian_W0 N ω) ≤ bulkEdge c + ε}) atTop (𝓝 1)

theorem singleTableLaw_of_moments [∀ N, IsProbabilityMeasure (μ N)] {c : ℝ} (hc : 0 < c)
    (m : SpikedModel μ n d) (hreg : m.Regime c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.GeneralNoise ν) :
    m.SingleTableLaw c
```

The three read one law `ν` and no Gaussian input. Under the law `noiseMatrix ν n d` itself,
with no model wrapper, the same statement is `Edge.tendsto_measure_opNorm_edge`, which the
target above transports through `(hG N).measure_eq`:

```lean
theorem tendsto_measure_opNorm_edge {c : ℝ} (hc : 0 < c) {ν : Measure ℝ} (hν : NoiseLaw ν)
    {n d : ℕ → ℕ} (hn : ∀ N, 0 < n N) (hd : ∀ N, 0 < d N)
    (hdtop : Tendsto d atTop atTop) (hcN : Tendsto (fun N => (n N : ℝ) / d N) atTop (𝓝 c))
    {a : ℝ} (ha : 0 < a) (ha8 : a < 1 / 8) :
    ∀ ε > 0, Tendsto (fun N => (noiseMatrix ν (n N) (d N))
      {Z | ‖Z‖ ^ 2 / (d N : ℝ) ≤ bulkEdge c + ε}) atTop (𝓝 1)
```

`a` is the truncation exponent: the proof truncates each entry at `d^(1/2 - a)` and the
target fixes `a = 1/16`. The one input of the chain, the Furedi-Komlos count now in
`RMT/General/Edge/CountBound.lean` (`cellCard_mul_le_pos`), was proved 2026-09-10 evening;
section 7.1 names it, and `docs/SORRIES.md` has no row, so the three theorems above and the
five `*_of_moments` forms are proved unconditionally. The chain is
`RMT/General/Edge/{Defs, Arith, Trunc, Sparse, Trace, Compare, Gaussian, Count, Code, Dyck,
CountBound, Excess, Markov, Sup}`, fourteen files, and the plan they follow is route C of
`notes/stage3_edge.md` (the campaign log is `notes/STAGE3_CAMPAIGN.md`).

**Hypotheses in words.**

- `[∀ N, IsProbabilityMeasure (μ N)]`: each `μ N` is a probability measure. Implicit in the
  paper.
- `hc : 0 < c`: the limiting aspect ratio is positive. This is `eq:RMT_limit`'s
  `c_i ∈ (0, ∞)`.
- `m : SpikedModel μ n d`: the model of `eq:model_rank1`.
- `hreg : m.Regime c`: `n N → ∞`, `d N → ∞` and `n N / d N → c`. This is `eq:RMT_limit`.
- `hG : m.GaussianNoise`: the unscaled noise matrix `Z N` has the law of a standard Gaussian
  matrix (`Defs.lean:168`), and the model's `E N` is `(√d)⁻¹ • Z N` (`Defs.lean:160`). The paper's `assum:general_noise` is wider; see modeling choice 1.

The general-noise form keeps `[∀ N, IsProbabilityMeasure (μ N)]`, `hc`, `m` and `hreg`. It
replaces `hG : m.GaussianNoise` by four binders.

- `hν : NoiseLaw ν`: the fixed law `ν` on `ℝ` is a probability measure with mean `0`, second
  moment `1` and a finite fourth moment (`Prob/NoiseLaw.lean`, quoted in section 2.6). One law
  serves every `N`.
- `hac : ν ≪ (volume : Measure ℝ)`: `ν` has a Lebesgue density. The paper does not ask for
  this; see below.
- `hG : m.GeneralNoise ν`: the unscaled noise matrix `Z N` has i.i.d. entries of law `ν` at
  every `N` (`Prob/NoiseLaw.lean:65`). At `ν = gaussianReal 0 1` it is `m.GaussianNoise`, by
  `Iff.rfl` (`SpikedModel.gaussianNoise_iff`).
- `hedge`: for every `ε > 0` the probability that the top eigenvalue of the downdated Gram
  matrix `W₀` is at most `bulkEdge c + ε` tends to `1`. This is the upper edge of the bulk,
  at the sharp constant `(1 + √c)²`.

**Hypotheses beyond the paper.** `[∀ N, IsProbabilityMeasure (μ N)]` is technical. The paper
works on an unnamed probability space, and the instance only fixes `μ N Set.univ = 1`, which
the events of probability tending to one need. `hc : 0 < c` is paper-implicit: `eq:RMT_limit`
(`main_paper.tex:262`) writes `n_i/d → c_i ∈ (0, ∞)`. `hreg : m.Regime c` is paper-implicit for
the same reason; it is `assum:main` (`main_paper.tex:271`) written as one limit statement.
`m : SpikedModel μ n d` is data, and it is the paper's `eq:model_rank1`. `hG : m.GaussianNoise`
is restrictive. `assum:general_noise` (`main_paper.tex:246`) asks only for i.i.d. entries with
mean `0`, unit variance after the `√d` scaling, and a fourth moment below a constant. The
Gaussian law is one member of that class, so what is lost is the four-moment universality.
Layer 1 takes the law as a hypothesis, so a four-moment proof replaces this one theorem and no
statement downstream. `singleTableLaw_of_general` (2026-09-09) is that proof, up to the edge.

**Hypotheses beyond the paper, general noise.** `hν : NoiseLaw ν` and `hG : m.GeneralNoise ν`
together are `assum:general_noise` (`main_paper.tex:246`) at one fixed law: i.i.d. entries,
mean `0`, unit variance after the `√d` scaling, and a fourth moment that is finite rather than
bounded by a constant. That is what Stage 1 buys, and it is the four-moment universality that
`hG : m.GaussianNoise` gives up. Two binders go beyond the paper. `hedge` is the upper edge of
the bulk. It is not paper-implicit in the sense of the binders above: the paper's proof takes
it from the literature, and the Gaussian chain proves it here by Sudakov-Fernique
(`tendsto_measure_lamMax_le`, `RMT/R3.lean:635`), which reads the Gaussian law. Stage 3 of
`notes/NONGAUSSIAN_SCOPE.md` named Bai-Yin as the route that would discharge it at a general
law, and `lamMax_W0_edge_of_general` (2026-09-10) is that discharge, by the moment method
with a truncation. `singleTableLaw_of_general` keeps `hedge` as a binder, so it stays an
honest implication; `singleTableLaw_of_moments` supplies the binder and keeps `hac` alone
beyond the paper.
`hac : ν ≪ (volume : Measure ℝ)` is restrictive: `assum:general_noise` allows an atomic law,
and a density excludes Rademacher entries. It pays for two steps. The first is the `topSimple`
field, which is false for Rademacher at `n = d = 2` and `θ = 0`, where the two eigenvalues of
`Xᵀ X` are `2 ± |b|` with `b = z₁₁z₁₂ + z₂₁z₂₂` and so coincide with probability `1/2`
(`notes/archive/prop_single_table_general.md` section 4). The second is the subcritical
`align` step, which reads the same field (`RMT/General/DelocAlign.lean:181`). The other three
fields hold without `hac`. Follow-up item F33 of `notes/FOLLOWUP_LIST.md` records the
Rademacher extension: weaken the field to "simple with probability tending to one", discharge
it by interlacing plus `hedge`, and drop `hac`; 1200 to 1800 Layer 1 lines.

**Modeling choices** (`notes/archive/prop_single_table.md`).

1. Gaussian noise in place of the paper's four-moment condition. This is user decision Q3 of
   `notes/FLAGGED.md`. The Layer 1 form is noise-free, so a four-moment discharge plugs in
   with no change to any downstream statement. `singleTableLaw_of_general` (Stage 1,
   2026-09-09) is that discharge, given the edge hypothesis `hedge` and a law with a density.
   It confirms the design: no Layer 1 statement changed, and the four representative
   corollaries of `General/Layer1.lean` (sections 3.2, 3.4, 3.5 and 5.4) are the Gaussian ones
   with the discharge swapped.
2. The overlap is the projector `‖P_top w‖²`, not `⟨v̂, w⟩²`. The two agree when the top
   eigenvalue is simple, which `topSimple` gives almost surely. This removes the sign
   ambiguity of `v̂` at no cost.
3. Part 2 is stated uniformly over unit `w ⊥ v` rather than for one deterministic sequence.
   Uniform implies sequential (`SingleTableLaw.deloc_seq`) and is what the cross-table lemma
   needs. It is proved from right-rotation invariance of the Gaussian model (`notes/FLAGGED.md`, D6).
4. `topSimple` is a field. It is stronger than anything the paper states, and it is proved
   for Gaussian noise by a resultant argument (`RMT/Simplicity.lean`), so it costs nothing
   here.
5. The threshold is `c < θ⁴` rather than `θ ≥ c^{1/4}`. Both branches give `0` at equality,
   and the equivalence uses the model field `0 ≤ θ` (`SpikedModel.hθ`).

**Proof sketch.** The proof is one case split on `c < θ⁴` inside `singleTableLaw_of_gaussian`
(`RMT/Full.lean:49`). Above the threshold the route is milestone L2a,
`singleTableLaw_of_gaussian_supercritical` (`RMT/Sup.lean:213`), reached through the tail shift
`singleTableLaw_of_gaussian_supercritical'` (`RMT/TailShift.lean:201`) that drops an auxiliary
`2 ≤ n N`. Item R0 splits one row off the noise and writes `Xᵀ X = W₀ + q qᵀ` with
`q = θ v + g` (`RMT/R0.lean:252`), so the top eigenvalue is the largest root of a secular
equation (`RMT/R4.lean:630`). The analytic core R5 is model-free: it consumes an edge bound and
six resolvent limits, bundled as `ResolventLimits`, and returns `align_tendstoInProb`
(`RMT/R5.lean:730`) and `lamMax_tendstoInProb` (`RMT/R5.lean:706`). Its route is an
intermediate value argument that localizes the secular root, then monotone brackets on the
quadratic forms. The two probabilistic inputs are proved here for the Gaussian law: the upper
edge `tendsto_measure_lamMax_le` (`RMT/R3.lean:635`), from Sudakov-Fernique on a net of the
unit spheres and one-sided Gaussian concentration, and the complex Marchenko-Pastur limits of
`RMT/R1.lean`, from a Stein identity (`Prob/GaussianAdapters.lean:195`) and a variance bound,
carried to the real axis by item T. Below the threshold the same interface, built by
`resolventLimits_of_gaussian'` (`RMT/TailShift.lean:190`), feeds
`align_tendstoInProb_of_subcritical` (`RMT/R6.lean:879`) and
`lamMax_tendstoInProb_of_subcritical'` (`RMT/TailShift.lean:260`). The last two fields hold in
both regimes: `delocUniform_of_gaussian` (`RMT/Symmetry.lean:594`) follows from right rotation
invariance of the Gaussian law, and `singleTableLaw_topSimple_of_gaussian`
(`RMT/Simplicity.lean:370`) is almost sure simplicity from a polynomial null set. This theorem
is itself the Layer 2 discharge, and every Gaussian facade in sections 3.2 to 3.7 calls it.
Sudakov-Fernique, Borell-TIS and the Gaussian concentration inequality are vendored in
`Vendor/COLT83/`, and Weyl and Courant-Fischer come from StatsMLlib; every random matrix step
in `RMT/` is proved in this project.

**Proof sketch, general noise.** The route is
`notes/archive/prop_single_table_general.md` sections 3 and 5, and it shares the model-free
analytic core with the Gaussian chain. Item R0 is replaced, not repaired: at a general law no
rotation makes the leave-one-out vector `g = Eᵀ u` independent of the remaining block, so
Stage 1 keeps the deterministic split `Xᵀ X = W₀ + q qᵀ` (`RMT/R0.lean`) and rewrites every
`g` form through Sherman-Morrison on the full noise Gram matrix `W = Eᵀ E`, which does have
i.i.d. rows, plus the companion resolvent. The six complex resolvent forms, bundled as the
structure `ResolventFormsC` (`RMT/General/Defs.lean`), then become rational expressions in
three isotropic forms. `resolventFormsC_of_general'` (`RMT/General/FormsGeneral.lean:148`)
proves them at four moments from two inputs: the isotropic law of `RMT/General/Iso.lean` and
`RMT/General/IsoMixed.lean`, a Chebyshev bound whose rate is uniform in the two directions
(`tendsto_isoRate`, `RMT/General/Iso.lean:2447`), built on the fourth moment of a linear form
(`Prob/LinFormMoments.lean`); and the trace law of `RMT/General/Trace.lean`, the
Marchenko-Pastur limit of the normalized trace at four moments. Every estimate runs at
`Im z > 0`, where `‖G(z)‖ ≤ 1/η` holds for every realization, so no estimate needs the edge;
item T (`RMT/T.lean`, model-free) buys the real axis. `hedge` enters once, as the `edge` field
of the seven-field `ResolventLimits` interface, in `resolventLimits_of_general`
(`RMT/General/Sup.lean:45`). The analytic core `RMT/R5.lean` is then used unchanged, so
`align` and `lamMax` above the threshold come from the same two theorems as in the Gaussian
proof. Below the threshold `lamMax` comes from `RMT/General/R3minus.lean:511` and `align` from
`align_tendstoInProb_of_subcritical_general'` (`RMT/General/DelocAlign.lean:122`). `topSimple`
is `singleTableLaw_topSimple_of_general` (`RMT/General/Simplicity.lean:128`), the polynomial
null set of the Gaussian item S at `noiseMatrix ν n d ≪ volume`, which is where `hac` is used.
`delocUniform` is unit G7, and it is not a port of a Gaussian argument. Its core is the window
bound `delocUniform_of_lamMax` (`RMT/General/DelocUniform.lean:203`), which says that on the
event `|λ_max - x| ≤ η` the overlap of the top eigenvector with a unit `w ⊥ v` is at most
`2 η Im (wᵀ G(x + iη) w)`, and the bridge identities of `RMT/General/FormsBridge.lean` reduce
`wᵀ G w` to forms whose bad events have `w`-free rates, which is what a supremum over `w`
needs. The two-regime wrapper is `delocUniform_of_general'`
(`RMT/General/DelocUniformSup.lean:107`). The window bound itself is unit D1a
(`RMT/General/Deloc.lean`, deterministic), and it serves two fields: `delocUniform` through
unit G7 and the subcritical `align` through unit D1c (`RMT/General/DelocAlign.lean`). Together
they replace `RMT/Symmetry.lean` (right-rotation invariance) and `RMT/R6.lean` (a conditional
rotation), neither of which survives a non-Gaussian law. No `hn2 : ∀ N, 2 ≤ n N` appears, because no row is split off, so
`RMT/TailShift.lean` is not used. By design
(`notes/archive/prop_single_table_general.md` choice 8) the general files re-prove the trace
law and the isotropic law from moment bounds instead of reusing `RMT/R1.lean` and
`RMT/R2.lean`, whose proofs rest on a Stein identity.

**Scope.** Convergence in probability, not almost surely. `lamMax` is a field of the structure
although the paper states it only inside the proof of `thm:theta_est`.
`singleTableLaw_of_gaussian` is unconditional. `singleTableLaw_of_general` is wider than it in
the noise law, which is the paper's own class, and narrower in two binders, `hedge` and `hac`,
so neither theorem implies the other. Both regimes are covered in both forms. The two chains
share the files `RMT/{MP7, R0, R4C, R5, ResolvDeriv, T}`, and of `RMT/R0.lean` the general
chain reads the deterministic content only, never its Gaussian block law.

**Audit trail.** `notes/archive/audit_scope_2026-08-29.md`,
`notes/archive/audit_numeric_2026-08-29.md` (closed forms against `theory_pred.R`),
`notes/archive/audit_l2_chain_2026-08-30.md`, `notes/archive/audit_mechanical_2026-08-31.md`.
Three external reviews examined the chain: `notes/archive/external_audit_2026-08-30.md`,
`external_audit_v2_review_2026-08-31.md`, `notes/archive/external_audit_v3_review_2026-09-01.md`.
The general-noise chain has its own trail: the plan note
`notes/archive/prop_single_table_general.md` (status `user OK` 2026-09-09, with the hypothesis
necessity scan of workflow rule 5 in its section 4, which drops each of the nine binders in
turn), the 21 unit reports `notes/archive/agent_reports/stage1_*.md`, and the numeric audit
`notes/archive/audit_nongaussian_forms_2026-09-09.md`, which recomputes the resolvent forms
and their rates at three laws with printed seeds (seed base 20260909).

### 3.2 `lem:delocalization`

**Paper.**

```latex
\begin{lemma} \label{lem:delocalization}
Let $\hat{v}_i$ be the top right singular vector of $X_i= \theta_i u_i v^\top + E_i$, $X_i \in \R^{n_i \times d}$ for $i=1,2$, where $X_1$ and $X_2$ are independent. Then, under \Cref{assum:general_noise,assum:main},
    \begin{equation*}
    |\langle \hat{v}_1, \hat{v}_2 \rangle|^2 \pto \beta_1^2 \beta_2^2.
    \end{equation*}
\end{lemma}
```

**Lean** (`SVDStack/Gram.lean`):

```lean
theorem lem_delocalization (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise)
    {i j : Fin M} (hij : i ≠ j) :
    TendstoInProb μ (fun N ω => ⟪m.vhat i N ω, m.vhat j N ω⟫_ℝ)
      (beta (m.tbl i).θ (c i) * beta (m.tbl j).θ (c j))
```

**Lean, general noise** (`General/Layer1.lean`, Stage 1 of `notes/NONGAUSSIAN_SCOPE.md`,
2026-09-09):

```lean
theorem lem_delocalization_of_general [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.JointGeneralNoise ν)
    (hedge : ∀ (i : Fin M), ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax ((m.tbl i).W0 N ω) ((m.tbl i).isHermitian_W0 N ω)
        ≤ bulkEdge (c i) + ε}) atTop (𝓝 1))
    {i j : Fin M} (hij : i ≠ j) :
    TendstoInProb μ (fun N ω => ⟪m.vhat i N ω, m.vhat j N ω⟫_ℝ)
      (beta (m.tbl i).θ (c i) * beta (m.tbl j).θ (c j))
```

The form without `hedge`, `lem_delocalization_of_moments` (`General/Layer1.lean`,
2026-09-10), has the same binders minus the edge; it is `lem_delocalization_of_general`
applied to `lamMax_W0_edge_of_general` at every table:

```lean
theorem lem_delocalization_of_moments [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.JointGeneralNoise ν)
    {i j : Fin M} (hij : i ≠ j) :
    TendstoInProb μ (fun N ω => ⟪m.vhat i N ω, m.vhat j N ω⟫_ℝ)
      (beta (m.tbl i).θ (c i) * beta (m.tbl j).θ (c j))
```

**Hypotheses in words.**

- `m : MultiTableModel μ M n d`: `M` tables of `eq:model_rank1` with one shared `v`.
- `c : Fin M → ℝ`: the limiting aspect ratios. The Layer 1 form takes them as data; `law`
  ties them to the tables.
- `law : ∀ i, (m.tbl i).SingleTableLaw (c i)`: the single-table law of every table.
- `hI : m.IndepNoise`: the tables are independent, each with some probability law on its
  noise matrix (`SVDStack/Gram.lean`; `JointGaussianNoise.indepNoise` gives it for Gaussian
  tables). The proof needs cross-table independence for a Fubini step over the product
  measure and nothing else about the noise; see modeling choice 2.
- `{i j : Fin M}` and `hij : i ≠ j`: the two tables compared, which must differ.

**Hypotheses beyond the paper.** None of the five binders of the Layer 1 form adds a condition
the paper lacks. `m` and `c` are data: the `M` tables of `eq:model_rank1` and their aspect
ratios. `law : ∀ i, (m.tbl i).SingleTableLaw (c i)` is paper-implicit, because it is
`prop:single_table` itself, which the paper applies to each table under `assum:general_noise`
and `assum:main`. `hI : m.IndepNoise` is stated in the paper's own lemma ("where `X_1` and
`X_2` are independent", `main_paper.tex:1112`), and it also follows from the entries of
`assum:general_noise` being i.i.d. across the table index (`main_paper.tex:246`). `hij : i ≠ j`
is the paper's `i = 1, 2`. The Layer 1 form asks for no `0 < c_i` and for no probability
measure, so it is weaker in its binders than every other statement of this section. The
finite-`d` companion below adds three: `hG : m.JointGaussianNoise` is restrictive in the same
way as in 3.1, `hd : 2 ≤ d N` is technical (at `d N = 1` the constant divides by zero and the
event is not small), and `hε : 0 < ε` is technical.

**Hypotheses in words, general noise.** `lem_delocalization_of_general` discharges `law`
itself, so it takes the model hypotheses `hc` and `hreg` in its place, and it replaces
`hI : m.IndepNoise` by the three noise binders `hν : NoiseLaw ν`,
`hac : ν ≪ (volume : Measure ℝ)` and `hG : m.JointGeneralNoise ν`: the `M` noise matrices are
independent, each with i.i.d. entries of the one fixed law `ν`, which has mean `0`, second
moment `1`, a finite fourth moment and a Lebesgue density (section 2.6 for the definitions,
section 3.1 for the two binders beyond the paper). The edge is taken **per table**: `hedge`
here is one hypothesis for each of the `M` tables, at that table's own aspect ratio `c i`.
The proof is one line: `m.lem_delocalization` fed by `SpikedModel.singleTableLaw_of_general`
per table and by `JointGeneralNoise.indepNoise` (section 3.4). No Layer 1 statement changes.

**Modeling choices** (`notes/archive/thm_svd_stack_general.md`, D1 of `notes/FLAGGED.md`).

1. The Lean statement is the **signed** limit `⟪v̂_i, v̂_j⟫ → β_i β_j`, not the squared one.
   The sign convention lives in `m.vhat`, which fixes `⟪v̂_i, v⟫ ≥ 0`. The squared form
   follows.
2. The noise predicate is `IndepNoise`, which says what the paper says: "independent"
   (`notes/archive/L4_indepnoise.md`). The Gaussian facades below keep `hG` and pass
   `hG.indepNoise`.

**Proof sketch.** The proof follows the paper's split
`v̂_iᵀ v̂_j = v̂_iᵀ v vᵀ v̂_j + v̂_iᵀ (I - v vᵀ) v̂_j` (`main_paper.tex:1123`) in
`lem_delocalization` (`SVDStack/Gram.lean:221`). The first term is a product of two signed
alignments. `align_inner` (`SVDStack/Gram.lean:200`) turns the `align` field of
`SingleTableLaw`, a limit of the squared overlap, into `⟪v̂_i, v⟫ → β_i`, by a square root on
the almost sure event where the top eigenvalue is simple (`overlap_eq_inner_sq`,
`Spectral.lean:516`). The second term is bounded, not computed. `abs_inner_perpOf_le`
(`SVDStack/DelocDir.lean:218`) replaces the component of `v̂_j` orthogonal to `v` by its normalized
direction, and `overlap_ge_inner_sq` (`Spectral.lean:528`) bounds the inner product by the
overlap of table `i` at that direction. That direction is random and reads table `j`, so the
uniform `delocUniform` field is applied through `measure_deloc_le` (`SVDStack/Gram.lean:187`).
Its proof `measure_deloc_le_of_pi_indep` (`SVDStack/DelocDir.lean:387`) is one Fubini step over the
product law, and this is the only place `hI` is used. The two pieces are combined by the
closure lemmas of `Prob/TendstoInProb.lean`, in particular `TendstoInProb.of_le`
(`Prob/TendstoInProb.lean:247`). Layer 2 discharges `law` with `singleTableLaw_of_gaussian`
(`RMT/Full.lean:49`) and `hI` with `JointGaussianNoise.indepNoise` (`SVDStack/Gram.lean:181`),
in the facades of `SVDStack/Main.lean`. Nothing here comes from StatsMLlib; the Mathlib input
is the product measure and Fubini theory.

**Scope.** Same as the paper. The general-noise form is the paper's noise class with a
Lebesgue density and one edge hypothesis per table added, as in section 3.1.

**Audit trail.** `notes/archive/audit_scope_2026-08-29.md`; the second half of
`notes/archive/thm_svd_stack_general.md`; the general form,
`notes/archive/prop_single_table_general.md` section 2.

**Beyond the paper (F1).** A finite-`d` tail bound, with an explicit constant, for
the cross term `⟪v̂_i, (I − vvᵀ) v̂_j⟫` (the second term of the paper's split at
`main_paper.tex:1123`), at the one
place the paper leans on an external asymptotic theorem (`main_paper.tex:1138`, Theorem 1 part
2 of `liu2023asymptotic`).

**Lean** (`SVDStack/DelocFinite.lean`):

```lean
theorem measure_inner_perpOf_ge_le (m : MultiTableModel μ M n d) (hG : m.JointGaussianNoise)
    {i j : Fin M} (hij : i ≠ j) (N : ℕ) (hd : 2 ≤ d N) {ε : ℝ} (hε : 0 < ε) :
    μ N {ω | ε ≤ |⟪m.vhat i N ω, perpOf ((m.tbl i).v N) (m.vhat j N ω)⟫_ℝ|}
      ≤ ENNReal.ofReal (1 / (ε ^ 2 * ((d N : ℝ) - 1)))

theorem measure_inner_perpOf_ge_le_toReal (m : MultiTableModel μ M n d)
    (hG : m.JointGaussianNoise) {i j : Fin M} (hij : i ≠ j) (N : ℕ) (hd : 2 ≤ d N) {ε : ℝ}
    (hε : 0 < ε) :
    (μ N {ω | ε ≤ |⟪m.vhat i N ω, perpOf ((m.tbl i).v N) (m.vhat j N ω)⟫_ℝ|}).toReal
      ≤ 1 / (ε ^ 2 * ((d N : ℝ) - 1))
```

**Hypotheses in words.** `hG : m.JointGaussianNoise`: independent Gaussian tables, stronger
than the `IndepNoise` the qualitative limit above needs; right-rotation invariance of the
Gaussian law is what gives the explicit constant. `hd : 2 ≤ d N`: at `d N = 1` the bound
divides by zero, and the event is not small (every unit vector is `± v`). `hε : 0 < ε`. The
bound reads none of `θ`, `n_i`, `c` or `M`, which is the point of the result.

**Scope.** New result beyond the paper: no paper statement gives a rate here, and the
positioning of Section 1 (`main_paper.tex:143`, `:150`) against non-asymptotic analyses costs
nothing to extend. See `notes/archive/F_batch_2026-09-05.md` for the proof plan (three finite-`N` facts
already in the tree, chained at `η = ε²`) and `notes/FOLLOWUP_LIST.md` item F1.

### 3.3 `lem:entrywise_conv_eigenvec`

**Paper.**

```latex
\begin{lemma} \label{lem:entrywise_conv_eigenvec}
    Suppose that $A_\beta$ has a unique largest eigenvalue $\lambda_1$ with $x =v_{\text{max}}(A_\beta)\in \R^M$, where without loss of generality $\langle x, \hat{x} \rangle \geq 0$. Then $\hat{x} = v_{\text{max}}(\tilde{V}\tilde{V}^\top)$ satisfies $(\hat{x})_m \pto x_m$ for all $m \in [M]$.
\end{lemma}
```

**Lean** (`SVDStack/Deterministic.lean`):

```lean
theorem lem_entrywise_conv_eigenvec (β : Fin M → ℝ)
    (G : (N : ℕ) → Ω N → Matrix (Fin M) (Fin M) ℝ)
    (hsymm : ∀ N ω, (G N ω).IsHermitian)
    (hconv : ∀ i j, TendstoInProb μ (fun N ω => G N ω i j) (Abeta β i j))
    (hsimple : TopSimple (Abeta β) (isHermitian_Abeta β))
    (x : EuclideanSpace ℝ (Fin M)) :
    TendstoInProb μ (fun N ω => ‖topProj (G N ω) (hsymm N ω) x‖ ^ 2)
      (⟪vMax (Abeta β) (isHermitian_Abeta β), x⟫_ℝ ^ 2)
```

The paper's signed coordinate form (`SVDStack/EntrywiseSigned.lean`; item L2 of the
external statement audit):

```lean
noncomputable def signPos (t : ℝ) : ℝ := if 0 ≤ t then 1 else -1

theorem lem_entrywise_conv_eigenvec_signed (β : Fin M → ℝ)
    (G : (N : ℕ) → Ω N → Matrix (Fin M) (Fin M) ℝ)
    (hsymm : ∀ N ω, (G N ω).IsHermitian)
    (hconv : ∀ i j, TendstoInProb μ (fun N ω => G N ω i j) (Abeta β i j))
    (hsimple : TopSimple (Abeta β) (isHermitian_Abeta β)) (m : Fin M) :
    TendstoInProb μ
      (fun N ω => signPos ⟪vMax (G N ω) (hsymm N ω), vMax (Abeta β) (isHermitian_Abeta β)⟫_ℝ
        * vMax (G N ω) (hsymm N ω) m)
      (vMax (Abeta β) (isHermitian_Abeta β) m)
```

**Hypotheses in words.** `G` is any sequence of symmetric random `M × M` matrices; `hconv` is
entrywise convergence in probability to `A_β`; `hsimple` is the paper's "unique largest
eigenvalue"; `x` is a fixed test vector. In the signed form, `m` is the coordinate and
`signPos ⟪x̂_N, x⟫ ∈ {1, -1}` is the sign that makes `⟪s_N x̂_N, x⟫ = |⟪x̂_N, x⟫| ≥ 0`, the
paper's "without loss of generality `⟨x, x̂⟩ ≥ 0`".

**Hypotheses beyond the paper.** `G` and `hsymm` are technical: the paper writes the one matrix
`Ṽ Ṽᵀ`, and Lean names an arbitrary symmetric random family so that the lemma is model-free.
Symmetry holds by construction for `Ṽ Ṽᵀ`, so the binder costs the consumer one
`isHermitian_mul_transpose_self`. `hconv` is paper-implicit: it is `lem:delocalization` read
entrywise, which step (iv) of the paper states at `main_paper.tex:1188`. `hsimple` is the
paper's "unique largest eigenvalue `λ_1`". `x` is technical: the first
conclusion is stated at a test vector, and the paper's coordinate form is the case
`x = v_max(A_β)`. In the signed form `m : Fin M` is the paper's `m ∈ [M]`, and the factor
`signPos ⟪x̂_N, x⟫` in the conclusion is technical, a written-out version of the paper's
"without loss of generality `⟨x, x̂⟩ ≥ 0`". No binder narrows the paper's claim.

**Modeling choices** (`notes/archive/thm_svd_stack_general.md`).

1. The conclusion of the first theorem is the projector overlap, not the entrywise limit "up
   to sign". The paper fixes the sign by fiat (`⟨x, x̂⟩ ≥ 0`); the projector form squares the
   sign away and is what `thm:svd_stack_general` consumes. The second theorem is the paper's
   coordinate form, with the sign written out as `signPos ⟪x̂_N, x⟫` (`1` when the inner
   product is `0`); it is a corollary of the first and has no consumer
   (`notes/archive/L2_entrywise_signed.md`).
2. The lemma is stated for an arbitrary `G`, not only for `Ṽ Ṽᵀ`. It is a special case of a
   model-free lemma `topProj_overlap_tendsto` in the same file.

**Proof sketch.** `lem_entrywise_conv_eigenvec` (`SVDStack/Deterministic.lean:360`) is the case
`A = A_β` of the model-free `topProj_overlap_tendsto` (`SVDStack/Deterministic.lean:295`). The
deterministic step is Davis-Kahan in projector form, `topProj_perturb`
(`LinAlg/TopProjPerturb.lean:243`): a gap `γ` for `A` and `‖B - A‖ ≤ δ` with `2 δ < γ` give
simplicity of `B` and `‖topProj B - topProj A‖ ≤ 4 δ / γ`. From it,
`norm_topProj_sq_continuousAt_of_topSimple` (`LinAlg/TopProjPerturb.lean:629`) states that
`‖topProj B x‖²` is continuous in the entries of `B` at a matrix with a simple top eigenvalue,
and `l2_opNorm_le_of_entries` (`LinAlg/TopProjPerturb.lean:365`) is the bridge from entries to
the operator norm. The probabilistic step is one application of the continuous mapping lemma
`TendstoInProbPi.comp_continuous` (`Prob/TendstoInProb.lean:383`) to the joint entrywise
convergence `tendstoInProbPi_entries` (`SVDStack/Deterministic.lean:273`). The gap that
Davis-Kahan needs comes from `abeta_gap` (`SVDStack/Deterministic.lean:179`), proved here by a
Rayleigh quotient at `e_i + e_j` against the second sorted eigenvalue. The signed form
`lem_entrywise_conv_eigenvec_signed` (`SVDStack/EntrywiseSigned.lean:96`) is a corollary.
`topSimple_whp_of_topSimple` (`SVDStack/EntrywiseSigned.lean:76`), through
`topSimple_whp_of_tendsto` (`SVDStack/Deterministic.lean:321`), makes the top eigenvalue of
`G_N` simple with probability tending to one; on that event the projector norm equals
`⟪x̂_N, x⟫²` (`norm_topProj_sq_eq_inner_sq`, `LinAlg/TopProjPerturb.lean:218`), so
`|⟪x̂_N, x⟫| →_p 1` and `‖s_N x̂_N - x‖² = 2 - 2 |⟪x̂_N, x⟫| →_p 0`. There is no Layer 2 step,
because the lemma reads no noise law; the consumer `thm_svd_stack_general` supplies `hconv`
from `gramEntries`. Weyl and Davis-Kahan enter through
`StatsMLlib.LinearAlgebra.Matrix.Perturbation`, and the projector forms of
`LinAlg/TopProjPerturb.lean` are derived from them here.

**Scope.** The paper's display, in both forms. The first theorem concludes
`‖P_top(G_N) x‖² →_p ⟨v_max(A_β), x⟩²` for a fixed test vector `x`, which is sign-free; the
second concludes `(s_N x̂_N)_m →_p x_m` for every coordinate `m`, with `x = v_max(A_β)` and
the sign `s_N` of `⟨x̂_N, x⟩`, which is the paper's statement. The second is derived from the
first: `topSimple_whp_of_topSimple` (same file; `topSimple_whp_of_tendsto` for `M ≥ 2`, every
`ω` for `M = 1`) makes the top eigenvalue of `G_N` simple with probability tending to `1`, on
that event `‖P_top(G_N) x‖² = ⟨x̂_N, x⟩²`, so `|⟨x̂_N, x⟩| →_p 1`,
`‖s_N x̂_N − x‖² = 2 − 2|⟨x̂_N, x⟩| →_p 0`, and so does every coordinate. Model-free
otherwise: the hypothesis is entrywise convergence in probability of `G_N`, with
no noise law. The route is Davis-Kahan, whose gap hypothesis is supplied by `abeta_gap`
(`λ_1(A_β) - λ_2(A_β) ≥ β_i β_j`).

**Audit trail.** `notes/archive/thm_svd_stack_general.md`,
`notes/archive/audit_l2_chain_2026-08-30.md`, `notes/archive/L2_entrywise_signed.md`.

### 3.4 `prop:stacksvd_general`

**Paper.**

```latex
\begin{prop}\label{prop:stacksvd_general}
    Under \Cref{assum:general_noise,assum:main}, the performance of \stacksvd is:
    \begin{equation*}
        | \langle \hat{v}_\stacksvd, v \rangle |^2 \pto \begin{cases}
            \frac{\|\theta\|_2^4-\|c\|_1}{ \|\theta\|_2^2(\|\theta\|_2^2 + 1)} & \text{ if } \|\theta\|_2^4 > \|c\|_1, \\
                   0 & \text{otherwise.}
               \end{cases}
    \end{equation*}
\end{prop}
```

**Paper, as amended (E10), no change to the environment.** E10.1 rewrites the sentence at `main_paper.tex:403` that calls Theorem 2.3 of `10.3150/19-BEJ1129` a heteroscedastic result; the reference assumes white noise, so the sentence says instead that the argument is adapted, and names the three steps that change.

**Lean** (`StackSVD/Main.lean`, since F28, 2026-09-08), Layer 1 and Gaussian:

```lean
theorem prop_stacksvd_general [NeZero M] (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (law : m.stack.SingleTableLaw (∑ i, c i)) :
    TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N))
      (stackSVDLimit (fun i => (m.tbl i).θ) c)

theorem prop_stacksvd_general_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N))
      (stackSVDLimit (fun i => (m.tbl i).θ) c)
```

The paper's own inner-product form is `prop_stacksvd_general_inner` and
`prop_stacksvd_general_inner_gaussian`, which take any selection `vhat` of a unit top
eigenvector of the stack Gram matrix.

**Lean, general noise** (`General/Layer1.lean`, Stage 1 of `notes/NONGAUSSIAN_SCOPE.md`,
2026-09-09), the facade and the two bridges of `General/Layer1.lean`:

```lean
theorem prop_stacksvd_general_of_general [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.JointGeneralNoise ν)
    (hedge : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax (m.stack.W0 N ω) (m.stack.isHermitian_W0 N ω)
        ≤ bulkEdge (∑ i, c i) + ε}) atTop (𝓝 1)) :
    TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N))
      (stackSVDLimit (fun i => (m.tbl i).θ) c)

theorem stack_law_general [NeZero M] (m : MultiTableModel μ M n d) {ν : Measure ℝ}
    [SigmaFinite ν] (h : m.JointGeneralNoise ν) : m.stack.GeneralNoise ν

theorem JointGeneralNoise.indepNoise {m : MultiTableModel μ M n d} {ν : Measure ℝ}
    [IsProbabilityMeasure ν] (hG : m.JointGeneralNoise ν) : m.IndepNoise
```

The form without `hedge`, `prop_stacksvd_general_of_moments` (`General/Layer1.lean`,
2026-09-10), has the same binders minus the edge; it is `prop_stacksvd_general_of_general`
applied to `lamMax_W0_edge_of_general` at the stack:

```lean
theorem prop_stacksvd_general_of_moments [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.JointGeneralNoise ν) :
    TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N))
      (stackSVDLimit (fun i => (m.tbl i).θ) c)
```

**Hypotheses in words.**

- `m : MultiTableModel μ M n d`: the `M` tables.
- `c : Fin M → ℝ`: the limiting aspect ratios.
- `law : m.stack.SingleTableLaw (∑ i, c i)` (Layer 1 only): the single-table law of the
  **stack**, at the parameters `(‖θ‖₂, ‖c‖₁)`. This is `prop:single_table` read once, on one
  matrix; see modeling choice 1.
- `hc : ∀ i, 0 < c i` (Gaussian only): `eq:RMT_limit`.
- `hreg : ∀ i, (m.tbl i).Regime (c i)` (Gaussian only): the regime of every table.
- `hG : m.JointGaussianNoise` (Gaussian only): independent Gaussian tables.
- `[∀ N, IsProbabilityMeasure (μ N)]` (Gaussian only): each `μ N` is a probability measure.

**Hypotheses beyond the paper.** `[NeZero M]` is technical: the shared `v` is read off table
`0`, so `Fin M` must be nonempty for the statement to elaborate (`SVDStack/Defs.lean:480` for
the same pattern). `law : m.stack.SingleTableLaw (∑ i, c i)` is paper-implicit in the Layer 1
form, because it is `prop:single_table` at the parameters `(‖θ‖₂, ‖c‖₁)`, which is what the
paper's proof applies to the stacked matrix. `hc : ∀ i, 0 < c i` and
`hreg : ∀ i, (m.tbl i).Regime (c i)` are paper-implicit, from `eq:RMT_limit`
(`main_paper.tex:262`) and `assum:main` (`main_paper.tex:271`). `hG : m.JointGaussianNoise` has
two halves. The Gaussian half is restrictive against the four-moment class of
`assum:general_noise` (`main_paper.tex:246`), exactly as in 3.1;
`prop_stacksvd_general_of_general` removes it, at the price of a Lebesgue density and one edge
hypothesis at the stack. The independence half is
paper-implicit, since `assum:general_noise` makes the entries i.i.d. across the table index `i`
as well. `[∀ N, IsProbabilityMeasure (μ N)]` is technical.

**Hypotheses in words, general noise.** `hν : NoiseLaw ν`, `hac : ν ≪ (volume : Measure ℝ)`
and `hG : m.JointGeneralNoise ν` replace `hG : m.JointGaussianNoise` (section 2.6 for the
definitions, section 3.1 for the two binders beyond the paper). The edge is taken **at the
stack**, one hypothesis, at the aspect ratio `∑ i, c i`, and not per table. That is the same
reduction as modeling choice 1 below: `stack_law_general` shows that the vertical stack of `M`
independent blocks with i.i.d. entries of law `ν` is one matrix with i.i.d. entries of law `ν`,
so the general single-table law applies to `m.stack` directly. It is the twin of `stack_law`
with the same proof, a `Measure.pi_eq` computation that never reads the Gaussian law.
`JointGeneralNoise.indepNoise` is the second bridge, the twin of
`JointGaussianNoise.indepNoise`; sections 3.2 and 3.5 use it for their Fubini step.

**Modeling choices** (`notes/archive/prop_stacksvd_general.md`).

1. The Gaussian reduction. With unit weights and independent Gaussian tables the stack is
   itself one rank-one spiked matrix with `θ_stack = ‖θ‖₂` and aspect ratio `‖c‖₁`
   (`MultiTableModel.stack`, `stack_law`, `stack_regime`). So this result needs no second
   black box; it is `prop:single_table` at `(‖θ‖₂, ‖c‖₁)`. `stack_law_general` is the same
   reduction at a general law, which is why `prop_stacksvd_general_of_general` takes one edge
   hypothesis at the stack and none per table.
2. The Layer 1 form takes no `IsProbabilityMeasure` instance, so it also holds for the zero
   measure, where it says nothing. That is a property of an implication, not a gap; the
   Gaussian form pins the measure through `GaussianNoise`
   (`notes/archive/audit_mechanical_2026-08-31.md`, findings 1 to 3).

**Proof sketch.** The reduction does the work here, and it is exact rather than asymptotic.
`MultiTableModel.stack` (`StackSVD.lean:157`) builds one `SpikedModel` on the vertical
concatenation, with signal `θ_stack = ‖θ‖₂` (`stack_theta`, `StackSVD.lean:171`) and the same
right singular vector (`stack_v`, `StackSVD.lean:175`). Two transport lemmas carry the
hypotheses to it. `stack_regime` (`StackSVD.lean:217`) gives `(∑_i n_i)/d → ∑_i c_i`, and
`stack_law` (`StackSVD.lean:292`) shows that a vertical stack of independent Gaussian blocks
has the law of one Gaussian matrix with `∑_i n_i` rows. `prop_stacksvd_general`
(`StackSVD/Main.lean:77`) is then three lines: read the `align` field of the stacked
`SingleTableLaw`, rewrite `θ_stack`, and turn `betaSq (√T) C` into the paper's closed form in
`T = ‖θ‖₂²` and `C = ‖c‖₁` with `betaSq_sqrt` (`StackSVD.lean:54`). There is no second
probabilistic step and no perturbation argument. Layer 2 discharges the stacked law in
`prop_stacksvd_general_gaussian` (`StackSVD/Main.lean:93`) by one call to
`singleTableLaw_of_gaussian` (`RMT/Full.lean:49`) on `m.stack`. Beyond the measure theory, the
only Mathlib input is the product law bookkeeping of `Prob/GaussianMatrix.lean` that
`stack_law` uses.

**Scope.** Same as the paper, in probability. The general-noise form is the paper's noise
class with a Lebesgue density and one edge hypothesis at the stack added, as in section 3.1;
it needs no per-table edge, because the stack is again one matrix with i.i.d. entries.

**Audit trail.** `notes/archive/audit_scope_2026-08-29.md` (section 3.2 item 1 asked for the
Gaussian link, which `prop_stacksvd_general_gaussian` supplies),
`notes/archive/audit_numeric_2026-08-29.md`, `notes/archive/audit_mechanical_2026-08-31.md`.

### 3.5 `thm:svd_stack_general`

**Paper.**

```latex
\begin{prop}
\label{thm:svd_stack_general}
    Under \Cref{assum:general_noise,assum:main}, when $\beta_2 > 0$, \svdstack satisfies:
    \begin{equation*}
        |\langle \hat{v}_\svdstack, v \rangle|^2 \pto
        \frac{\left(\beta^\top v_{\text{max}}(A_\beta)\right)^2}{\lambda_\text{max}\left(A_\beta\right)}.
    \end{equation*}
    When $\beta_1 = 0$, this inner product converges in probability to 0.
\end{prop}
```

**Paper, as amended (E9), no change to the environment.** E9.1 adds the `M ≥ 2` qualifier to the sentence at `main_paper.tex:376` that calls `β_2 > 0` necessary and sufficient; at `M = 1` the matrix `A_β = [1]` has a unique largest eigenvalue and `β_2` does not exist.

**Lean** (`SVDStack/Main.lean`):

```lean
theorem thm_svd_stack_general [NeZero M] (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ i j, i ≠ j ∧ 0 < β i ∧ 0 < β j)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) (svdstackLimit β)

theorem thm_svd_stack_general_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ i j, i ≠ j ∧ 0 < β i ∧ 0 < β j)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) (svdstackLimit β)
```

Companions: `thm_svd_stack_general_inner` and `_inner_gaussian` (the paper's inner product),
`thm_svd_stack_general_zero` and `_zero_gaussian` (the `β_1 = 0` clause, limit `0`).

**Lean, general noise** (`General/Layer1.lean`, Stage 1 of `notes/NONGAUSSIAN_SCOPE.md`,
2026-09-09):

```lean
theorem thm_svd_stack_general_of_general [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ i j, i ≠ j ∧ 0 < β i ∧ 0 < β j)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.JointGeneralNoise ν)
    (hedge : ∀ (i : Fin M), ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax ((m.tbl i).W0 N ω) ((m.tbl i).isHermitian_W0 N ω)
        ≤ bulkEdge (c i) + ε}) atTop (𝓝 1)) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) (svdstackLimit β)
```

The form without `hedge`, `thm_svd_stack_general_of_moments` (`General/Layer1.lean`,
2026-09-10), has the same binders minus the edge; it is `thm_svd_stack_general_of_general`
applied to `lamMax_W0_edge_of_general` at every table:

```lean
theorem thm_svd_stack_general_of_moments [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ i j, i ≠ j ∧ 0 < β i ∧ 0 < β j)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.JointGeneralNoise ν) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) (svdstackLimit β)
```

**Hypotheses in words.**

- `hβdef` names `β_i = beta θ_i c_i`, so `β` is not free.
- `hthr : ∃ i j, i ≠ j ∧ 0 < β i ∧ 0 < β j` is the paper's `β_2 > 0`, written without an
  ordering. It gives the spectral gap `λ_1(A_β) - λ_2(A_β) ≥ β_i β_j > 0` through
  `abeta_gap`, and it forces `2 ≤ M`.
- `hI : m.IndepNoise` (Layer 1) supplies the cross-table independence that
  `lem_delocalization` needs, with arbitrary marginal laws; the Gaussian facade takes
  `hG : m.JointGaussianNoise` and passes `hG.indepNoise` (L4).

**Hypotheses beyond the paper.** `[NeZero M]` is technical and cannot be dropped from the
statement: `svdstackPerf` reads the shared `v` off table `0` (`SVDStack/Defs.lean:480`), so
`0 : Fin M` must elaborate. Semantically it adds nothing, since `hthr` already forces `2 ≤ M`
through `one_lt_card_fin_of_ne` (`SVDStack/Deterministic.lean:165`).
`hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)` is technical: it names `β`, so that the limit
`svdstackLimit β` reads as in the paper. `hc : ∀ i, 0 < c i` is paper-implicit (`eq:RMT_limit`,
`main_paper.tex:262`), and the proof uses it only to get `0 ≤ β_i < 1` through `beta_mem_Ico`
(`SVDStack/Deterministic.lean:60`). `hthr` is the paper's `β_2 > 0` without an ordering. `law`
is paper-implicit, being `prop:single_table` per table. `hI : m.IndepNoise` is paper-implicit
through the i.i.d. entries of `assum:general_noise` (`main_paper.tex:246`), and the paper
states it outright in `lem:delocalization`. In the Gaussian form `hreg` is paper-implicit, `hG`
is restrictive in its Gaussian half as in 3.1, and `[∀ N, IsProbabilityMeasure (μ N)]` is
technical. `thm_svd_stack_general_of_general` removes the Gaussian half, at the price of a
Lebesgue density and one edge hypothesis per table. The `_inner` companion adds `vhat`,
`hmem` and `hnorm`, all technical: they let the
reader name any unit top eigenvector, and no measurability of it is required.

**Hypotheses in words, general noise.** `hβdef`, `hthr` and `hc` are unchanged.
`hν : NoiseLaw ν`, `hac : ν ≪ (volume : Measure ℝ)` and `hG : m.JointGeneralNoise ν` replace
`hG : m.JointGaussianNoise`, and `hreg` stays because the theorem discharges `law` itself
(section 2.6 for the definitions, section 3.1 for the two binders beyond the paper). The edge
is taken **per table**, as in 3.2 and not as in 3.4: this result reads each table's own top
singular vector, so it needs the single-table law of each of the `M` tables, and `hedge` gives
one edge hypothesis for each at its own aspect ratio `c i`. The proof is one line:
`m.thm_svd_stack_general` fed by `SpikedModel.singleTableLaw_of_general` per table and by
`JointGeneralNoise.indepNoise`.

**Modeling choices** (`notes/archive/thm_svd_stack_general.md`).

1. `hthr` in unsorted form. The paper sorts `β`; sorting would need a permutation lemma on
   `A_β` and buys nothing.
2. `svdstackLimit` reads one arbitrary top eigenvector of `A_β` through `vMax`. When the top
   eigenvalue of `A_β` repeats, the value depends on the choice: at `β = (0.648, 0, 0)` the
   two readings are `0` and `0.419` (`notes/archive/audit_mechanical_2026-08-31.md`,
   finding 5). Every theorem that uses `svdstackLimit` carries simplicity, here from `hthr`.
3. `hc` replaces an earlier `0 ≤ β_i < 1` hypothesis, which follows from the formula for
   `betaSq`.
4. The simplicity the proof needs on the svdstack Gram matrix holds with probability tending
   to one (`topSimple_svdstackGram`), not almost surely. The per-table `topSimple` of
   `SingleTableLaw` is almost sure and is free, because item S proves it (`notes/FLAGGED.md`, D2).

**Proof sketch.** `thm_svd_stack_general` (`SVDStack/Main.lean:219`) follows the paper's four
steps. Step one is the Gram matrix: `gramEntries` (`SVDStack/Gram.lean:315`) gives entrywise
convergence in probability of `Ṽ Ṽᵀ` to `A_β`, with the off-diagonal from `lem_delocalization`
and the diagonal exactly `1`. Step two is the gap: `abeta_gap`
(`SVDStack/Deterministic.lean:179`) gives `λ_1(A_β) - λ_2(A_β) ≥ β_i β_j > 0`, so
`topSimple_of_gap` (`LinAlg/TopProjPerturb.lean:187`) makes the top eigenvalue of the limit
simple. Step three works with the paper's closed form `(xᵀ Ṽ v)² / λ_max(Ṽ Ṽᵀ)` at
`x = v_max(Ṽ Ṽᵀ)` (`svdstackPerfClosed`, `SVDStack/Defs.lean:221`), and splits it into a
denominator and a numerator. The denominator converges by Weyl, in the form
`lamMax_tendstoInProb` (`SVDStack/Deterministic.lean:280`). The numerator converges at the
fixed vector `β` by `lem_entrywise_conv_eigenvec` (`SVDStack/Deterministic.lean:360`), and the
random vector `Ṽ v` is exchanged for `β` by the deterministic bound
`abs_norm_sq_topProj_sub_le` (`SVDStack/Main.lean:66`) together with `align`
(`SVDStack/Gram.lean:434`). Step four transfers the quotient back to `svdstackPerf`: the two
agree on a good event (`svdstackPerf_eq_of_goodEvent`, `SVDStack/Main.lean:171`) whose
complement has vanishing measure (`tendsto_measure_not_goodEvent`, `SVDStack/Gram.lean:329`),
and `TendstoInProb.of_tendsto_measure_ne_of_tendsto` (`Prob/TendstoInProb.lean:277`) closes it.
The paper's inner-product companion uses the same device with `topSimple_svdstackGram`
(`SVDStack/Gram.lean:412`). Layer 2 discharges `law` and `hI` in
`thm_svd_stack_general_gaussian` (`SVDStack/Main.lean:382`), by `singleTableLaw_of_gaussian`
(`RMT/Full.lean:49`) per table and `JointGaussianNoise.indepNoise` (`SVDStack/Gram.lean:181`).
Weyl and Davis-Kahan come from StatsMLlib; the Rayleigh quotient handling, the good event and
every probabilistic step are proved here.

**Scope.** Same as the paper. The zero clause is a separate theorem, as in the paper's last
sentence. The general-noise form is the paper's noise class with a Lebesgue density and one
edge hypothesis per table added, as in section 3.1; the zero clause has no general-noise twin,
because `General/Layer1.lean` carries the four representative facades only.

**Audit trail.** `notes/archive/thm_svd_stack_general.md`,
`notes/archive/audit_scope_2026-08-29.md`, `notes/archive/audit_mechanical_2026-08-31.md`,
`notes/archive/audit_independent_D33_2026-09-02.md` (the `Rayleigh`/`Simple` facades).

### 3.6 `thm:simple_thm1` (cor. 1)

**Paper.**

```latex
\begin{cor}\label{thm:simple_thm1}
    Suppose $c_i=c_0$ and $\theta_i=\theta_0$ for all $i$.
    Under \Cref{assum:general_noise,assum:main}, the performance of \stacksvd and \svdstack is given by:
\begin{align*}
    | \langle \hat{v}_\svdstack, v \rangle |^2 &\pto \begin{cases}
              1 - \frac{c_0 +  \theta_0^2}{M \theta_0^4 + \theta_0^2 - (M-1)c_0} & \text{ if } \theta_0 \geq c_0^{1/4}, \\
              0 & \text{ otherwise.}
             \end{cases}\\
    | \langle \hat{v}_\stacksvd, v \rangle |^2 &\pto \begin{cases}
                1 - \frac{c_0 +  \theta_0^2}{M \theta_0^4 + \theta_0^2} & \text{ if } \theta_0 \geq M^{-1/4} c_0^{1/4}, \\
                0 & \text{ otherwise.}
            \end{cases}
\end{align*}
\end{cor}
```

**Paper, as amended (E5), no change to the environment.** E5.1 adds one sentence after the corollary, at `main_paper.tex:338`: the statement holds for every `M ≥ 1`, at `M = 1` the SVDstack branch is `prop:single_table`, and the two branches of each display agree at the threshold.

**Lean.** stacksvd half (`StackSVD/Main.lean`, since F28), svdstack half in two forms
(`SVDStack/Main.lean`, `SVDStack/Simple.lean`):

```lean
theorem thm_simple_thm1_stacksvd_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (θ₀ c₀ : ℝ) (hc : 0 < c₀)
    (hθ : ∀ i, (m.tbl i).θ = θ₀) (hreg : ∀ i, (m.tbl i).Regime c₀)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N))
      (if c₀ < (M : ℝ) * θ₀ ^ 4 then
          1 - (c₀ + θ₀ ^ 2) / ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2) else 0)

theorem thm_simple_thm1_svdstack_gaussian_full [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (θ₀ c₀ : ℝ) (hc : 0 < c₀)
    (hθ : ∀ i, (m.tbl i).θ = θ₀) (hreg : ∀ i, (m.tbl i).Regime c₀) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω)
      (if c₀ < θ₀ ^ 4 then
          1 - (c₀ + θ₀ ^ 2) / ((M : ℝ) * θ₀ ^ 4 + θ₀ ^ 2 - ((M : ℝ) - 1) * c₀) else 0)
```

`thm_simple_thm1_svdstack_gaussian` is the `2 ≤ M`, above-threshold case;
`thm_simple_thm1_svdstack_gaussian_full` covers every `M ≥ 1` and both branches. The scalar
identities are `Scalars.simple_thm1_stacksvd` and `Scalars.simple_thm1_svdstack`.

**Hypotheses in words.** `hθ` and `hreg` are the paper's "`c_i = c_0` and `θ_i = θ_0` for all
`i`". `hc : 0 < c₀` is `eq:RMT_limit`.

**Hypotheses beyond the paper.** `[NeZero M]` is technical, again because the shared `v` is
table `0`'s, and `[∀ N, IsProbabilityMeasure (μ N)]` is technical. `hc : 0 < c₀` is
paper-implicit, from `eq:RMT_limit` (`main_paper.tex:262`). `hθ : ∀ i, (m.tbl i).θ = θ₀` and
`hreg : ∀ i, (m.tbl i).Regime c₀` are the corollary's own "`c_i = c_0` and `θ_i = θ_0` for all
`i`", with `hreg` also carrying `assum:main`. `hG : m.JointGaussianNoise` is restrictive in its
Gaussian half, as in 3.1, and paper-implicit in its independence half. The `hM2 : 2 ≤ M` and
`hdet : c₀ < θ₀⁴` binders of the intermediate `thm_simple_thm1_svdstack_gaussian`
(`SVDStack/Main.lean:428`) do not appear in the `_full` form quoted above, which covers every
`M ≥ 1` and both branches.

**Modeling choices** (`notes/archive/agent_reports/layer1_leftovers.md`, choice 1).

1. The `if` form carries both branches in one statement, so the reader sees the threshold.
2. At `M = 1` the svdstack formula reduces to `β²`, consistent with `prop:single_table`; the
   proof of `_full` makes that case explicit. This is `notes/paper_edits.md` E5 item 1.

**Proof sketch.** Both halves are corollaries, so the work is scalar algebra and a case split,
not new probability. The stacksvd half `thm_simple_thm1_stacksvd_gaussian`
(`StackSVD/Main.lean:128`) instantiates `prop_stacksvd_general_gaussian`
(`StackSVD/Main.lean:65`) at the constant `c ≡ c₀` and `θ ≡ θ₀`, then rewrites the general
limit into the paper's closed form with `Scalars.simple_thm1_stacksvd` (`Scalars.lean:802`)
directly (F32, 2026-09-08, drops the local restatement `stackSVDLimit_const_eq`). The svdstack half splits on the threshold
inside `thm_simple_thm1_svdstack_gaussian_full` (`SVDStack/Simple.lean:91`). Above the
threshold and with `2 ≤ M` it calls `thm_simple_thm1_svdstack_gaussian`
(`SVDStack/Main.lean:428`), which instantiates `thm_svd_stack_general_gaussian`
(`SVDStack/Main.lean:382`) at the constant `β` and rewrites by `Scalars.simple_thm1_svdstack`
(`Scalars.lean:830`); the two distinct indices that `hthr` needs come from `2 ≤ M`, and `β > 0`
comes from `c₀ < θ₀⁴`. At `M = 1` the svdstack Gram matrix has one entry, so `svdstackPerf`
equals the single-table overlap (`svdstackPerf_eq_overlap_one`, `SVDStack/Simple.lean:56`) and
the closed form collapses to `betaSq θ₀ c₀` by one field identity. Below the threshold every
`β_i` is `0`, and the limit `0` comes from `thm_svd_stack_general_zero_gaussian`
(`SVDStack/Main.lean:398`). Layer 2 needs no new step: both halves inherit the Gaussian
discharge of 3.4 and 3.5, which rests on `singleTableLaw_of_gaussian` (`RMT/Full.lean:49`). The
scalar identities in `Scalars.lean` are closed by `field_simp` and `nlinarith`, with no Mathlib
input beyond ordered field lemmas.

**Scope.** The paper's thresholds `θ₀ ≥ c₀^{1/4}` and `θ₀ ≥ M^{-1/4} c₀^{1/4}` are the Lean
guards `c₀ < θ₀⁴` and `c₀ < M θ₀⁴`; both branches give `0` at equality, and the equivalence
uses `0 ≤ θ₀` (the model field `SpikedModel.hθ`).

**Audit trail.** `notes/archive/agent_reports/layer1_leftovers.md` (choice 1),
`notes/archive/audit_independent_D33_2026-09-02.md` section 2 (verdict: pass with findings, 0
must-fix). The label `cor.1` versus `cor.2` was fixed as question Q6 of `notes/FLAGGED.md`.

### 3.7 `cor.2` (binary weighting)

**Paper.**

```latex
\begin{cor}\label{cor.2}
    Under \Cref{assum:general_noise,assum:main} the performance of binary-weighted \stacksvd, where tables with $\theta_i^4 \le c_i$ are discarded, is given by (assuming at least one table is above the threshold of detectability):
    \begin{align*}
        % \CS := \{ i : \theta_i^4 \geq c \} \\
        | \langle \hat{v}_\stacksvd, v \rangle |^2 \pto \frac{\left(\sum_{i \in \CS} \theta_i^2\right)^2 - \sum_{i \in \CS} c_i}{\left(\sum_{i \in \CS} \theta_i^2\right)^2 + \sum_{i \in \CS} \theta_i^2},\quad\text{where}\quad \CS := \{ i : \theta_i^4 > c_i \}.
    \end{align*}
\end{cor}
```

**Paper, as amended (E9), no change to the environment.** E9.2 rewrites the subset optimization at `main_paper.tex:439` to `442`, which sits after the corollary: the maximum runs over nonempty subsets and each subset carries its own indicator, so the `0/0` maximand at `S = ∅` and the negative value at a general subset both go away (R1 and R7).

The paper then optimizes over all subsets.

**Lean** (`StackSVD/Weighted.lean`, since F28):

```lean
theorem stackPerfW_binary_tendsto_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (S : Finset (Fin M)) [NeZero S.card] (c : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω)
      (Scalars.binaryStackSVDLimit S (fun i => (m.tbl i).θ) c)

theorem exists_binary_tendsto_max_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    ∃ S : Finset (Fin M), S.Nonempty ∧
      TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω)
        (Scalars.binaryStackSVDLimitMax (fun i => (m.tbl i).θ) c)
```

The scalar half is `Scalars.binaryStackSVDLimit_eq` and
`binaryStackSVDLimit_eq_stackSVDLimit`.

**Hypotheses in words.** `S` is any nonempty subset, not only the detectable set. The paper
says the same in the sentence after the corollary.

**Hypotheses beyond the paper.** `[NeZero M]` is technical, for the shared `v` of table `0`,
and `[∀ N, IsProbabilityMeasure (μ N)]` is technical. `[NeZero S.card]` is the paper's own
"assuming at least one table is above the threshold of detectability", stated here for a
general `S`. `hc : ∀ i, 0 < c i` and `hreg : ∀ i, (m.tbl i).Regime (c i)` are paper-implicit,
from `eq:RMT_limit` (`main_paper.tex:262`) and `assum:main` (`main_paper.tex:271`).
`hG : m.JointGaussianNoise` is restrictive in its Gaussian half and paper-implicit in its
independence half, as in 3.4. `S : Finset (Fin M)` is data, and it widens the paper's claim
rather than narrowing it: the paper fixes `S` to the detectable set, and the Lean statement
holds at every nonempty subset. The existence form `exists_binary_tendsto_max_gaussian` carries
the same binders without `S` and `[NeZero S.card]`.

**Modeling choices** (`notes/archive/class_A.md`, `notes/archive/agent_reports/class_A_scalars.md`).

1. The result is stated for an arbitrary nonempty `S`, and the paper's detectable set is one
   instance (used in section 5.1). The subset maximum is a separate existence statement,
   because `binaryStackSVDLimitMax` is a maximum over a finite family. That maximum
   (`Scalars.lean:1112`) is a `sup'` over every subset, the empty one included (value `0`),
   as in the paper's `2^{[M]}`; `Scalars.exists_nonempty_eq_max` (`Scalars.lean:1173`, which
   needs `[NeZero M]`) shows that a nonempty subset attains it.
2. The model half is `MultiTableModel.restrict S`, a genuine sub-collection, with
   `restrict_jointGaussianNoise` transporting the Gaussian law. So the corollary is
   `prop:stacksvd_general` on the sub-collection, not a new limit law.

**Proof sketch.** The corollary is `prop:stacksvd_general` on a sub-collection, and the only
new work is to show that the sub-collection is again a model with the same noise law.
`MultiTableModel.restrict` (`StackSVDWeighted.lean:1303`) builds the model on the `S.card` kept
tables. `restrict_jointGaussianNoise` (`StackSVDWeighted.lean:1417`) transports the joint
Gaussian law across two reindexings, a `Finset` restriction of a product measure and an
equivalence of index types, with the measure-preserving maps of Mathlib.
`heteroLaw_binary_of_gaussian` (`StackSVDWeighted.lean:1431`) then reads the stacked
single-table law of the restricted model, exactly as in 3.4, and
`heteroLaw_binary_of_singleTableLaw` (`StackSVDWeighted.lean:1370`) repackages it as a
`HeteroLaw` (`StackSVDWeighted.lean:1268`) at the binary weights.
`stackPerfW_binary_tendsto_gaussian` (`StackSVD/Weighted.lean:187`) is then two lines: take the
`align` field, and rewrite the weighted scalar limit at binary weights with `Scalars.Lw_binary`
(`StackSVDWeighted.lean:627`). The paper's closed form in `S` is
`Scalars.binaryStackSVDLimit_eq` (`Scalars.lean:174`), and
`binaryStackSVDLimit_eq_stackSVDLimit` (`Scalars.lean:147`) identifies it with
`prop:stacksvd_general` read on `S`. The subset maximum is a separate statement because the
maximum is attained, not computed: `Scalars.exists_nonempty_eq_max` (`Scalars.lean:1173`) picks
a nonempty maximizer over the finite family of subsets, and
`exists_binary_tendsto_max_gaussian` (`StackSVD/Weighted.lean:199`) applies the previous
theorem at it. Layer 2 discharges the law by `singleTableLaw_of_gaussian` (`RMT/Full.lean:49`)
on the restricted stack, through `stack_regime` (`StackSVD.lean:217`) and `stack_law`
(`StackSVD.lean:292`). No perturbation theory enters, so nothing from StatsMLlib is used in
this section.

**Scope.** Same as the paper. The strict detectable set uses `c_i < θ_i⁴`, so a table exactly
on the threshold is discarded; this matters for `remark:svd_outperform_stack`
(`notes/paper_edits.md` E6).

**Audit trail.** `notes/archive/class_A.md`, `notes/archive/audit_weighted_2026-08-30.md`,
`notes/archive/audit_global_2026-08-31.md`.

## 4. Weighted results

### 4.1 `thm:stacksvd_weighted`

**Paper.**

```latex
\begin{thm}
    \label{thm:stacksvd_weighted} Under \Cref{assum:general_noise,assum:main}, \stacksvd is optimally weighted as
    \begin{equation*}
    w_i\opt \propto \frac{\theta_i}{\sqrt{\theta_i^2 + c_i}},
\end{equation*}
which yields performance
\begin{equation*}
    (v^\top \hat{v}_\stacksvd)^2 \pto \gamma\opt \quad \text{ the unique solution $x \in (0,1)$ of } \sum_{i=1}^M \theta_i^4 \frac{1-x}{c_i + x\theta_i^{2}} = 1,
\end{equation*}
as long as $\sum_i \theta_i^4/ c_i > 1$, otherwise $\gamma\opt = 0$.
\end{thm}
```

**Paper, as amended (E8).** The theorem gains the hypothesis that at least one `θ_i` is positive. At `θ ≡ 0` the stated weights are all zero, the weighted stack is the zero matrix, and the claim fails; the exact replacement is E8.1 in `notes/paper_edits.md`.

```latex
\begin{thm}
    \label{thm:stacksvd_weighted} Under \Cref{assum:general_noise,assum:main}, and with $\theta_i > 0$ for at least one $i$, \stacksvd is optimally weighted as
    \begin{equation*}
    w_i\opt \propto \frac{\theta_i}{\sqrt{\theta_i^2 + c_i}},
\end{equation*}
which yields performance
\begin{equation*}
    (v^\top \hat{v}_\stacksvd)^2 \pto \gamma\opt \quad \text{ the unique solution $x \in (0,1)$ of } \sum_{i=1}^M \theta_i^4 \frac{1-x}{c_i + x\theta_i^{2}} = 1,
\end{equation*}
as long as $\sum_i \theta_i^4/ c_i > 1$, otherwise $\gamma\opt = 0$.
\end{thm}
```

**Paper, as amended (E9 and E10), no change to the environment.** E9.5 repairs the inequality at `main_paper.tex:1612`, which points the wrong way, and states instead the sign equivalence that the next item uses; E10.2 names the adaptation of Theorem 2.3 of `10.3150/19-BEJ1129` from white noise to the block variance profile at `main_paper.tex:1372`.

**Lean.** Layer 1 (`StackSVD/Weighted.lean`, since F28) and Gaussian (`RMT/Het/Sup.lean`):

```lean
theorem thm_stacksvd_weighted [NeZero M] (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hθ : ∃ i, (m.tbl i).θ ≠ 0)
    (law : m.HeteroLaw (Scalars.optWstack (fun i => (m.tbl i).θ) c) c) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c)

theorem thm_stacksvd_weighted_gaussian [NeZero M]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hθ : ∃ i, (m.tbl i).θ ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c)
```

The optimality half is bundled in:

```lean
theorem thm_stacksvd_weighted_gaussian_opt [NeZero M]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hθ : ∃ i, (m.tbl i).θ ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c)
    ∧ ∀ w : Fin M → ℝ, (∃ i, w i ≠ 0) →
        TendstoInProb μ (fun N ω => m.stackPerfW w N ω)
          (Scalars.Lw (fun i => (m.tbl i).θ) c w)
        ∧ Scalars.Lw (fun i => (m.tbl i).θ) c w
            ≤ Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c
```

`thm_stacksvd_weighted_inner` is the paper's `(vᵀ v̂)²`; `thm_stacksvd_weighted_general` is
the limit at any weight vector. `thm_stacksvd_weighted_gaussian_margin`, `_margin_inner`,
`_margin_smul` are the interim theorems under the Sudakov-Fernique margin, kept as a record.

**Hypotheses in words.**

- `hc : ∀ i, 0 < c i`: `eq:RMT_limit`.
- `hθ : ∃ i, θ_i ≠ 0`: at `θ = 0` the optimal weight vector is zero, the weighted stack is
  the zero matrix, and the claim is false. Not in the paper; in the paper as amended (E8).
  The Layer 1 proof does not use it (`law.align` is rewritten with `Scalars.L_optW_eq hc`,
  `StackSVDWeighted.lean:1003`); it works in the Gaussian form, where
  `Scalars.optWstack_ne_zero hc hθ` supplies the nonzero weight vector that
  `heteroLaw_of_gaussian` needs.
- `law`: the heteroscedastic law at the optimal weights.
- `hreg`, `hG`: the regime and the joint Gaussian law.

**Hypotheses beyond the paper.** `[NeZero M]` is technical: `stackPerfW` reads the shared `v`
off table `0` (`StackSVDWeighted.lean:1251`), so `Fin M` must be nonempty for the statement to
elaborate. `hc : ∀ i, 0 < c i` is paper-implicit, from `eq:RMT_limit` (`main_paper.tex:262`).
`hθ : ∃ i, θ_i ≠ 0` is restrictive (finding E8): it removes the one parameter point `θ = 0`, where the
paper's own weights `w_i⋆ ∝ θ_i/√(θ_i² + c_i)` are the zero vector, the weighted stack is the
zero matrix, and the `v_max` of `stacksvd.def` (`main_paper.tex:301`) is an arbitrary unit
vector, so the paper's `γ_opt = 0` fails there. Every `θ ≠ 0` stays covered, in both regimes.
`law : m.HeteroLaw w⋆ c` (Layer 1 only) is paper-implicit: the `align` field is
`eq:weighted_norm` (`main_paper.tex:1416`) at `b = v`, reduced to `L(w)`, and the `topSimple`
field is the simple top eigenvalue that the definition of `v̂_stacksvd` presumes.
`hreg : ∀ i, (m.tbl i).Regime (c i)` (Gaussian only) is paper-implicit, `assum:main`
(`main_paper.tex:271`). `hG : m.JointGaussianNoise` (Gaussian only) has two halves, as in 3.4.
The Gaussian half is restrictive against the four-moment class of `assum:general_noise`
(`main_paper.tex:246`); it costs every non-Gaussian law of that class (user decision Q3). The
independence half is paper-implicit, since `assum:general_noise` makes the entries i.i.d.
across the table index `i` as well. No `IsProbabilityMeasure` instance appears in either form:
at Layer 1 nothing pins the measure (modeling choice 2), and in the Gaussian form `hG` gives
it.

**Modeling choices** (`notes/archive/thm_stacksvd_weighted.md`).

1. The limit is a **total** function `Scalars.stackSVDLimitW`, with `0` below the threshold
   `∑ θ_i⁴/c_i > 1`. The paper's two cases are inside it. `Lw θ c 0 = 0` is a convention.
2. `HeteroLaw.align` is stated at `L(w)` (`main_paper.tex:1434`), one deterministic step past
   the BEJ1129 display it cites (`:1414`). The Gaussian proof proves that field directly, so
   the bridge is part of the proof, not a restatement of a black box (`notes/FLAGGED.md`, D21).
3. The optimality claim `Lw θ c w ≤ stackSVDLimitW θ c` is a separate scalar theorem
   `Scalars.L_le_opt`, proved by a simplex optimization with Cauchy-Schwarz. The bundle
   `_opt` states both halves at once.

**Proof sketch.** `MultiTableModel.stackW` (`StackSVDWeighted.lean:1151`) builds the weighted
stack as one `SpikedModel`, and `stackW_X` (`StackSVDWeighted.lean:1183`) checks that block `i`
of it is `w_i X_i`. Layer 1 is then two lines: `thm_stacksvd_weighted`
(`StackSVD/Weighted.lean:53`) reads the `align` field of `HeteroLaw`
(`StackSVDWeighted.lean:1268`) and rewrites `Scalars.Lw θ c w⋆` into `Scalars.stackSVDLimitW θ c`
with `Scalars.L_optW_eq` (`StackSVDWeighted.lean:1003`). That rewrite carries the paper's
appendix computation. `Scalars.assumption4_and_L_optW_of_thr` (`StackSVDWeighted.lean:864`) runs
the change of variables of `main_paper.tex:1595`: at the root `r` of
`∑_i θ_i⁴(1-r)/(c_i + rθ_i²) = 1` it sets `γ₁ = 1/(1-r)`, checks the secular equation through
`Scalars.sum_secular_eq_neg_one` (`StackSVDWeighted.lean:274`), and gets `L(w⋆) = r`.
`Scalars.assumption4_optW_iff` (`StackSVDWeighted.lean:1016`) identifies `eq:assumption4` at
`w⋆` with the threshold `1 < ∑_i θ_i⁴/c_i`, which is how the two branches of `stackSVDLimitW`
line up with the paper's two cases. The optimality half `Scalars.L_le_opt`
(`StackSVDWeighted.lean:845`) puts `p_j = θ_j²w_j²/(γ₁ - w_j²)` on the simplex, bounds the
constraint sum below by Cauchy-Schwarz in `Scalars.one_le_gW_Lw`
(`StackSVDWeighted.lean:742`), and closes with `Scalars.le_stackSVDLimitW` (`Scalars.lean:442`).
The paper's inner-product form `thm_stacksvd_weighted_inner` (`StackSVD/Weighted.lean:66`)
squeezes over `HeteroLaw.topSimple` through `overlap_eq_inner_sq`. Layer 2 discharges the whole
structure for Gaussian noise in `heteroLaw_of_gaussian` (`RMT/Het/Sup.lean:205`): `topSimple` is
`heteroLaw_topSimple_of_gaussian` (`RMT/Het/Simplicity.lean:270`), and `align` splits on
`eq:assumption4`, above it `align_tendstoInProb_het_of_assumption4` (`RMT/Het/R5het.lean:753`)
over the resolvent limits of `resolventLimitsHet_of_edge` (`RMT/Het/Sup.lean:146`), below it
`align_tendstoInProb_het_subcritical` (`RMT/Het/R6het.lean:396`). Both branches take the exact
upper edge from `heteroEdge_of_gaussian` (`RMT/Het/EdgeSharp.lean:785`), a sharp
Sudakov-Fernique bound at `γ⋆ = -1/s⋆`, so `thm_stacksvd_weighted_gaussian`
(`RMT/Het/Sup.lean:292`) is one line, with `Scalars.optWstack_ne_zero` (`RMT/Het/Sup.lean:88`)
for `w⋆ ≠ 0`. Only the Gaussian comparison and concentration inputs come from outside this
project: `sudakov_fernique` and Borell-TIS are the vendored Mathlib backport in
`Vendor/COLT83/Mathlib/Probability/`, while the resolvent chain of `RMT/Het/` is proved here.

**Scope.** Same as the paper, with `hθ` added. `heteroLaw_of_gaussian` covers both regimes
with no margin condition; the exact edge is `heteroEdge_of_gaussian`.

**Audit trail.** `notes/archive/thm_stacksvd_weighted.md`,
`notes/archive/audit_weighted_stacksvd_2026-08-30.md`,
`notes/archive/audit_het_skeleton_2026-08-31.md`,
`notes/archive/audit_het_chain_2026-09-01.md`, `notes/archive/heteroedge_sharp.md`. Question
Q11 and decisions D22, D27, D28 record the route change to the sharp edge.

### 4.2 `thm:svdstack_weighted`

**Paper.**

```latex
\begin{thm}
    \label{thm:svdstack_weighted}
    Under \Cref{assum:general_noise,assum:main}, an optimal weighting of \svdstack is
    \begin{equation}\label{svdstack.weight}
        w_i \opt = \theta_i \sqrt{\frac{\theta_i^2+1}{\theta_i^2+c_i}} \mathds{1} \left\{ \theta_i^4 > c_i\right\},
    \end{equation}
    which when $\beta_{1}>0$ yields performance
    \begin{equation*}
        |\langle v, \hat{v}_\svdstack \rangle|^2 \pto \frac{S}{S+1}, \quad \text{ where } \quad S = \sum_i \frac{\beta_i^2}{1-\beta_i^2}.
    \end{equation*}
\end{thm}
```

**Paper, as amended (E1).** The optimality claim becomes a bound over every nonzero weighting, deterministic or data dependent, with no eigengap condition, and the stated weighting attains it; the exact replacement is E1.1 in `notes/paper_edits.md`.

```latex
\begin{thm}
    \label{thm:svdstack_weighted}
    Under \Cref{assum:general_noise,assum:main}, set
    \begin{equation*}
        S = \sum_i \frac{\beta_i^2}{1-\beta_i^2}.
    \end{equation*}
    No weighting of \svdstack, deterministic or data-dependent, beats $S/(S+1)$: for every $\varepsilon > 0$,
    \begin{equation*}
        \mathbb{P}\left( \sup_{w \neq 0} |\langle v, \hat{v}_\svdstack(w) \rangle|^2 > \frac{S}{S+1} + \varepsilon \right) \to 0,
    \end{equation*}
    where the supremum runs over all nonzero $w \in \R^M$. The weighting
    \begin{equation}\label{svdstack.weight}
        w_i \opt = \theta_i \sqrt{\frac{\theta_i^2+1}{\theta_i^2+c_i}} \mathds{1} \left\{ \theta_i^4 > c_i\right\},
    \end{equation}
    attains this bound when $\beta_{1}>0$, and is then an optimal weighting:
    \begin{equation*}
        |\langle v, \hat{v}_\svdstack(w\opt) \rangle|^2 \pto \frac{S}{S+1}.
    \end{equation*}
    (When $\beta_1 = 0$, $S = 0$ and the bound already gives $|\langle v, \hat{v}_\svdstack(w) \rangle|^2 \pto 0$ for every nonzero $w$.)
\end{thm}
```

**Paper, as amended (E1 and E9), no change to the environment.** E1.3 puts the uniform bound into the proof at `main_paper.tex:1284`, before the per-weight limit; E1.2 adds the bound to Case 3 at `:1223`, where the paper gives none; E9.4 restricts the `w̃ ≤ 1` condition at `:1331` to the undetectable tables, since `w̃_i = 1/√(1-β_i²) > 1` on a detectable table, so no weighting satisfies both conditions (R3).

**Lean** (`SVDStack/Weighted.lean`):

```lean
theorem thm_svdstack_weighted [NeZero M] (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ k, 0 < β k)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerfW (optW β) N ω) (svdstackLimitOpt β)

theorem thm_svdstack_weighted_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ k, 0 < β k)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerfW (optW β) N ω) (svdstackLimitOpt β)

theorem thm_svdstack_weighted_paper [NeZero M] (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ k, 0 < β k)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hI : m.IndepNoise) :
    TendstoInProb μ
      (fun N ω => m.svdstackPerfW (fun i => paperW (m.tbl i).θ (c i)) N ω)
      (svdstackLimitOpt β)
```

`thm_svdstack_weighted_paper` is the literal `eq:svdstack.weight` weighting, where an
undetectable table gets weight `0`, not `optW β i = 1`; the limit is the same. The optimality
bundle is:

```lean
theorem thm_svdstack_weighted_gaussian_opt [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hthr : ∃ k, 0 < β k)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerfW (optW β) N ω) (svdstackLimitOpt β)
    ∧ ∀ w : Fin M → ℝ, (∃ k, w k ≠ 0) → (∀ i, 0 ≤ w i) →
        TopSimple (AbetaW w β) (isHermitian_AbetaW w β) →
        TendstoInProb μ (fun N ω => m.svdstackPerfW w N ω) (svdstackLimitW w β)
        ∧ svdstackLimitW w β ≤ svdstackLimitOpt β
```

`thm_svdstack_weighted_gaussian_opt_full` (`SVDStack/Rayleigh.lean`) adds a bound that is
uniform in `w`: with probability tending to one, no weight vector beats `S/(S+1) + ε`. Its
signature carries `hthr : ∃ k, 0 < β k` for the packaging with the attainment clause; the
uniform conjunct itself, `svdstackPerfW_uniform_bound` (`SVDStack/Rayleigh.lean:542`), needs
only `hc`, `hβdef`, `law` and `hI`. That
closes the paper's optimality claim without an eigengap condition, which is
`notes/paper_edits.md` E1.

**Hypotheses in words.**

- `m`, `c`, `β`: the tables, the aspect ratios, and the per-table overlaps.
- `hc : ∀ i, 0 < c i`: `eq:RMT_limit`.
- `hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)`: names `β`, so it is not free.
- `hthr : ∃ k, 0 < β k`: the paper's `β_1 > 0`. One detectable table is enough, unlike the
  unweighted `thm_svd_stack_general`, which needs two. There is no `2 ≤ M`: at `M = 1` the
  statement is `prop:single_table`.
- `law : ∀ i, (m.tbl i).SingleTableLaw (c i)` (Layer 1) or `hreg` and
  `[∀ N, IsProbabilityMeasure (μ N)]` (Gaussian): the random matrix theory input.
- `hI : m.IndepNoise` (Layer 1): cross-table independence with arbitrary marginals, which
  `lem_delocalization` needs; the Gaussian facades take `hG : m.JointGaussianNoise` and pass
  `hG.indepNoise`. See modeling choice 3.

**Hypotheses beyond the paper.** `[NeZero M]` is technical, for the reason of 3.4: table `0`
names the shared `v`. `[∀ N, IsProbabilityMeasure (μ N)]` (Gaussian only) is technical.
`hc : ∀ i, 0 < c i` is paper-implicit, `eq:RMT_limit` (`main_paper.tex:262`).
`hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)` is technical: it names `β` so that the statement can
quote the paper's `β_i`, and it constrains no model. `hthr : ∃ k, 0 < β k` is paper-implicit, the
theorem's own `β₁ > 0`. `law : ∀ i, (m.tbl i).SingleTableLaw (c i)` (Layer 1 only) is
paper-implicit: it is `prop:single_table` for each table, which the paper's proof cites at
`main_paper.tex:1284`. `hI : m.IndepNoise` (Layer 1 only) is paper-implicit, since
`assum:general_noise` (`main_paper.tex:246`) makes the entries i.i.d. across the table index; it
is the weakest form of that half, because the marginals stay arbitrary.
`hreg : ∀ i, (m.tbl i).Regime (c i)` (Gaussian only) is paper-implicit, `assum:main`
(`main_paper.tex:271`). `hG : m.JointGaussianNoise` (Gaussian only) splits as in 4.1: the
Gaussian half is restrictive against `assum:general_noise` and costs every non-Gaussian law of
that four-moment class (user decision Q3), the independence half is paper-implicit. Inside the
optimality quantifier of `thm_svdstack_weighted_gaussian_opt`, `∀ i, 0 ≤ w i` is restrictive
against the theorem's own `w ∈ R^M` (`main_paper.tex:453`; the estimator `svdstack.def` at
`:308` takes any `w`), while `main_paper.tex:300` writes `w ∈ R^M_{≥0}`; the uniform bound of
`_opt_full` has no sign condition. `TopSimple (AbetaW w β)` is paper-implicit,
the proviso "provided that `W A_β W` has a unique largest eigenvalue" of `main_paper.tex:1284`;
`thm_svdstack_weighted_gaussian_opt_full` removes both. The base theorem
`thm_svdstack_weighted_general` asks neither `0 ≤ w i` nor `2 ≤ M`, so on those two points the
Lean form is broader than the paper.

**Modeling choices** (`notes/archive/thm_svdstack_weighted.md`).

1. Two weightings are stated: the rewritten `optW β = 1/√(1-β_i²)` and the paper's literal
   `paperW`. `svdstackLimitW_paperW` identifies the limits.
2. The optimality quantifier ranges over nonnegative `w` with a simple top eigenvalue of
   `A_{β,w}`. The uniform Rayleigh bound of `SVDStack/Rayleigh.lean` removes both conditions
   at the cost of an `ε` and a "with probability tending to one".
3. `hI : m.IndepNoise` sits in the Layer 1 signature because `lem_delocalization` needs
   cross-table independence; nothing else about the noise enters at Layer 1 (L4).

**Proof sketch.** The estimator has two stages, so the proof reduces the `d × d` problem to an
`M × M` one, along the display at `main_paper.tex:1282`. On an almost sure event
`svdstackGramW_eq` (`SVDStack/Weighted.lean:544`) rewrites `∑_i w_i² P_i` as `Ṽ_wᵀ Ṽ_w`, and
`svdstackPerfW_eq_of_goodEventW` (`SVDStack/Weighted.lean:689`) turns the performance into
`‖P_top(Ṽ_w Ṽ_wᵀ) Ṽ_w v‖² / λ_max(Ṽ_w Ṽ_wᵀ)`. The probabilistic input is two limits: `align`
(`SVDStack/Gram.lean:434`) gives `Ṽ v → β` from `prop:single_table`, and `gramEntriesW`
(`SVDStack/Weighted.lean:563`) gives `Ṽ_w Ṽ_wᵀ → A_{β,w}` entrywise through `gramEntries`
(`SVDStack/Gram.lean:315`), whose off-diagonal is `lem_delocalization`
(`SVDStack/Gram.lean:221`) and is where `hI` enters. The deterministic step is perturbation at a
fixed size: `lamMax_tendstoInProb` (`SVDStack/Deterministic.lean:280`) is Weyl and
`topProj_overlap_tendsto` (`SVDStack/Deterministic.lean:295`) is Davis-Kahan, both from
StatsMLlib, pushed through the continuous mapping theorem of `Prob/TendstoInProb.lean`.
`thm_svdstack_weighted_general` (`SVDStack/Weighted.lean:902`) assembles the quotient, replaces
`Ṽ_w v` by its limit `Wβ` with the projector bound `abs_norm_sq_topProj_sub_le`, and transfers
off the good event by `tendsto_measure_not_goodEventW` (`SVDStack/Weighted.lean:599`). At the
optimal weights the eigengap proviso is free, because every `w_l²(1 - β_l²)` equals `1` there and
`abetaW_gap` (`SVDStack/Weighted.lean:118`) applies, so `thm_svdstack_weighted`
(`SVDStack/Weighted.lean:990`) needs one detectable table only. The value `S/(S+1)` comes from
the Sherman-Morrison step `abetaW_optW_eq` (`SVDStack/Weighted.lean:172`) and
`svdstackLimitW_eq_opt` (`SVDStack/Weighted.lean:203`), collected in `svdstackLimitW_optW`
(`SVDStack/Weighted.lean:349`); the `M = 1` case runs through
`thm_svdstack_weighted_of_card_le_one` (`SVDStack/Weighted.lean:883`), which is
`prop:single_table`. `thm_svdstack_weighted_paper` (`SVDStack/Weighted.lean:1036`) repeats this
at the paper's literal weights, identified by `svdstackLimitW_paperW`
(`SVDStack/Weighted.lean:386`), and `svdstackLimitW_le_opt` (`SVDStack/Weighted.lean:415`) is the
optimality, Cauchy-Schwarz in the `A_β` inner product. Layer 2 discharges the per-table law in
`thm_svdstack_weighted_gaussian` (`SVDStack/Weighted.lean:1337`), by
`SpikedModel.singleTableLaw_of_gaussian` (`RMT/Full.lean:49`) on each table and `hG.indepNoise`
for `hI`. The uniform bound is a separate argument: `svdstackPerfW_le_rowBound`
(`SVDStack/Rayleigh.lean:368`) bounds every nonzero weighting at once, on one event that does
not read `w`, by the squared norm of the projection of `v` on the row space of `Ṽ`;
`rowBound_eq_closed` (`SVDStack/Rayleigh.lean:416`) writes that as `gᵀ G⁻¹ g`, `rowBound_tendsto`
(`SVDStack/Rayleigh.lean:482`) sends it to `S/(S+1)`, and `svdstackPerfW_uniform_bound`
(`SVDStack/Rayleigh.lean:542`) is the union bound that `thm_svdstack_weighted_gaussian_opt_full`
(`SVDStack/Rayleigh.lean:590`) reports.

**Scope.** Matches the paper as amended (E1): the uniform bound is the amended statement.
Against the printed statement the uniform bound is stronger.

**Audit trail.** `notes/archive/thm_svdstack_weighted.md`,
`notes/archive/audit_weighted_2026-08-30.md` (decision D10),
`notes/archive/P1_rayleigh_bound.md`, `notes/archive/audit_independent_D33_2026-09-02.md` section 3
(verdict: pass with findings, 0 must-fix).

### 4.3 `lem:secular_equation`

**Paper.**

```latex
\begin{lemma}\label{lem:secular_equation}
    At least $n_i - 1$ eigenvalues of $R$ are equal to $w_i^2$ for $i \in [M]$, and the remaining eigenvalues are given by the roots of the secular equation
    \begin{equation*}
        f(\lambda) = 1 + \sum_{j=1}^M \frac{\theta_j^2 w_j^2}{w_j^2 - \lambda}.
    \end{equation*}
    Under condition \eqref{eq:assumption4}, the largest eigenvalue $\gamma_1$ is unique and a corresponding eigenvector $\xi_1$ is given by
    \begin{equation*}
        \xi_1 \propto (\Sigma - \gamma_1 I)^{-1} \tilde{u}_0.
    \end{equation*}
\end{lemma}
```

**Lean.** Two halves. The scalar half is in `StackSVDWeighted.lean`: `secular`,
`secular_strictMonoOn`, `IsGammaTop`, `sum_secular_eq_neg_one`. The spectral half is in
`Secular.lean`: `det_Rmat_sub` (the paper's determinant display),
`mem_spectrum_iff_secularDiag` (claim 2), `mulVec_eq_of_const_on` and
`card_sub_one_le_finrank_eigenspace` (claim 1), `root_above_unique`, `lamMax_Rmat_eq`,
`topSimple_Rmat`, `topSpace_Rmat_eq_span` (claim 3), plus the block forms
`secularDiag_stack`, `isGammaTop_iff_stack`, `lamMax_Rstack_eq_gammaTop`,
`topSpace_Rstack_eq_span`.

The spectral half, verbatim (`Secular.lean`; the section variables are `{nn : ℕ}`,
`{σ q : Fin nn → ℝ}` from line 288 on, and `{M : ℕ} {ν : Fin M → ℕ}` in the stack section;
`toOp` is the linear map of a matrix on `EuclideanSpace ℝ (Fin nn)`):

```lean
noncomputable def Rmat (σ q : Fin nn → ℝ) : Matrix (Fin nn) (Fin nn) ℝ :=
  Matrix.diagonal σ + Matrix.vecMulVec q q

noncomputable def secularDiag (σ q : Fin nn → ℝ) (lam : ℝ) : ℝ :=
  1 + ∑ i, q i ^ 2 / (σ i - lam)

theorem det_Rmat_sub (σ q : Fin nn → ℝ) {lam : ℝ} (hlam : ∀ i, σ i ≠ lam) :
    (Rmat σ q - lam • (1 : Matrix (Fin nn) (Fin nn) ℝ)).det
      = (∏ i, (σ i - lam)) * secularDiag σ q lam

theorem mem_spectrum_iff_secularDiag (σ q : Fin nn → ℝ) {lam : ℝ} (hlam : ∀ i, σ i ≠ lam) :
    lam ∈ spectrum ℝ (Rmat σ q) ↔ secularDiag σ q lam = 0

theorem card_sub_one_le_finrank_eigenspace (σ q : Fin nn → ℝ) {t : ℝ} {S : Finset (Fin nn)}
    (hS : ∀ i ∈ S, σ i = t) :
    S.card - 1 ≤ Module.finrank ℝ (Module.End.eigenspace (toOp (Rmat σ q)) t)

theorem root_above_unique {z₁ z₂ : ℝ}
    (h₁ : (⨆ i, σ i) < z₁) (h₂ : (⨆ i, σ i) < z₂)
    (e₁ : secularDiag σ q z₁ = 0) (e₂ : secularDiag σ q z₂ = 0) : z₁ = z₂

theorem lamMax_Rmat_eq {lam : ℝ} (hlam : (⨆ i, σ i) < lam)
    (e : secularDiag σ q lam = 0) :
    lamMax (Rmat σ q) (isHermitian_Rmat σ q) = lam

theorem topSimple_Rmat {lam : ℝ} (hlam : (⨆ i, σ i) < lam)
    (e : secularDiag σ q lam = 0) : TopSimple (Rmat σ q) (isHermitian_Rmat σ q)

theorem topSpace_Rmat_eq_span {lam : ℝ} (hlam : (⨆ i, σ i) < lam)
    (e : secularDiag σ q lam = 0) :
    topSpace (Rmat σ q) (isHermitian_Rmat σ q)
      = Submodule.span ℝ
          {(WithLp.toLp 2 fun i => q i / (σ i - lam) : EuclideanSpace ℝ (Fin nn))}

noncomputable def sigmaStack (ν : Fin M → ℕ) (w : Fin M → ℝ) : Fin (∑ i, ν i) → ℝ :=
  fun p => w (finSigmaFinEquiv.symm p).1 ^ 2

noncomputable def u0Stack (θ w : Fin M → ℝ) (u : (i : Fin M) → Fin (ν i) → ℝ) :
    Fin (∑ i, ν i) → ℝ :=
  fun p => θ (finSigmaFinEquiv.symm p).1 * w (finSigmaFinEquiv.symm p).1 *
    u (finSigmaFinEquiv.symm p).1 (finSigmaFinEquiv.symm p).2

theorem secularDiag_stack (θ w : Fin M → ℝ) {u : (i : Fin M) → Fin (ν i) → ℝ}
    (hu : ∀ i, ∑ k, u i k ^ 2 = 1) (lam : ℝ) :
    secularDiag (sigmaStack ν w) (u0Stack θ w u) lam = Scalars.secular θ w lam

theorem lamMax_Rstack_eq_gammaTop (θ w : Fin M → ℝ) {u : (i : Fin M) → Fin (ν i) → ℝ}
    (hu : ∀ i, ∑ k, u i k ^ 2 = 1) (hν : ∀ i, 0 < ν i) (hM : 0 < M)
    (hg : ∃ g, Scalars.IsGammaTop θ w g) :
    lamMax (Rmat (sigmaStack ν w) (u0Stack θ w u))
        (isHermitian_Rmat (sigmaStack ν w) (u0Stack θ w u)) = Scalars.gammaTop θ w
      ∧ TopSimple (Rmat (sigmaStack ν w) (u0Stack θ w u))
        (isHermitian_Rmat (sigmaStack ν w) (u0Stack θ w u))

theorem topSpace_Rstack_eq_span (θ w : Fin M → ℝ) {u : (i : Fin M) → Fin (ν i) → ℝ}
    (hu : ∀ i, ∑ k, u i k ^ 2 = 1) (hν : ∀ i, 0 < ν i) (hM : 0 < M)
    (hg : ∃ g, Scalars.IsGammaTop θ w g) :
    topSpace (Rmat (sigmaStack ν w) (u0Stack θ w u))
        (isHermitian_Rmat (sigmaStack ν w) (u0Stack θ w u))
      = Submodule.span ℝ
          {(WithLp.toLp 2 fun p => u0Stack θ w u p /
              (sigmaStack ν w p - Scalars.gammaTop θ w) :
            EuclideanSpace ℝ (Fin (∑ i, ν i)))}
```

How the block reads against the paper. `Rmat (sigmaStack ν w) (u0Stack θ w u)` is the
paper's `R = Σ + ũ₀ ũ₀ᵀ` with `Σ = diag(w_i² I_{n_i})` and block `i` of `ũ₀` equal to
`θ_i w_i u_i`; `secularDiag_stack` collapses the block sum to the paper's `f(λ)`. Claim 1
("at least `n_i − 1` eigenvalues equal `w_i²`") is `card_sub_one_le_finrank_eigenspace` at
`t = w_i²` and `S` the set of all stacked indices with `σ_p = t`: when two tables share a
weight the level set is the union of their blocks, so a coincidence of weights raises the
multiplicity bound to `∑ n_i − 1` over the tied blocks; the paper's per-table count is the
special case of disjoint level sets. Claim 2 ("the remaining eigenvalues are the roots") is
`mem_spectrum_iff_secularDiag` away from the diagonal values, with `det_Rmat_sub` the
determinant identity behind it. Claim 3 is `root_above_unique` (the root above `max σ` is
unique), `lamMax_Rmat_eq` (it is the top eigenvalue), `topSimple_Rmat` (it is simple) and
`topSpace_Rmat_eq_span` (the eigenvector is `(Σ − γ₁ I)⁻¹ ũ₀`, written entrywise as
`q_i/(σ_i − λ)`), all under `q ≠ 0` and the existence of a root above `sup σ`; the stack
forms `lamMax_Rstack_eq_gammaTop` and `topSpace_Rstack_eq_span` instantiate them at
`γ₁ = Scalars.gammaTop θ w` under `∃ g, IsGammaTop θ w g`, which is where `eq:assumption4`
enters. (External statement audit finding 15.)

**Hypotheses in words.** All deterministic, on a finite matrix at a fixed size; no model, no
measure, no limit. The scalar declarations take `θ w : Fin M → ℝ` and, where the top root must
exist, `∃ g, Scalars.IsGammaTop θ w g`. The spectral declarations take `σ q : Fin nn → ℝ` and
build the matrix `Rmat σ q` with its Hermitian proof `isHermitian_Rmat`. No declaration of
the section takes `Scalars.Assumption4` or a `c`; the paper's `eq:assumption4` enters the
downstream theorems (section 4.1) as `Scalars.Assumption4 θ c w`, which unfolds to
"the outlier root exists, and `∑_i c_i w_i⁴/(γ₁ - w_i²)² < 1`".

**Hypotheses beyond the paper.** The section is deterministic, so most binders are the paper's
own conditions written out. The abstract setting `{nn : ℕ}` with `{σ q : Fin nn → ℝ}` is more
general than the paper, not less: the diagonal is arbitrary, and the paper's block `Σ` is the
instance `sigmaStack`. `hlam : ∀ i, σ i ≠ lam` in `det_Rmat_sub` and
`mem_spectrum_iff_secularDiag` is paper-implicit, since the paper's `f(λ)` divides by
`w_j² - λ` (`main_paper.tex:1381`). `hS : ∀ i ∈ S, σ i = t` in
`card_sub_one_le_finrank_eigenspace` is paper-implicit, the block form of `Σ`
(`main_paper.tex:1360`). `hlam : (⨆ i, σ i) < lam` with `e : secularDiag σ q lam = 0`, and the
scalar `hg : ∃ g, Scalars.IsGammaTop θ w g`, are paper-implicit: they are `γ₁ > max_i w_i²` and
the secular equation of `main_paper.tex:1370`, the existence half of `eq:assumption4`.
`hu : ∀ i, ∑ k, u i k ^ 2 = 1` is paper-implicit, because block `i` of the paper's `ũ₀` is
`θ_i w_i u_i` with `‖u_i‖ = 1`. `hν : ∀ i, 0 < ν i` and `hM : 0 < M` are technical: an empty block
or an empty table index leaves no top eigenvalue to name, and the paper has `n_i ≥ 1` and
`M ≥ 1`. `hq : q ≠ 0` and `hnn : 0 < nn` in `root_above_unique`, `lamMax_Rmat_eq`,
`topSimple_Rmat` and `topSpace_Rmat_eq_span` left the four signatures on 2026-09-07 (F30). Both
binders are derivable: `secularDiag σ q lam` is the constant `1` when `q = 0`, so the root
hypothesis `e` already excludes `q = 0`, and a nonzero `q` on `Fin nn` forces `0 < nn`. The two
new lemmas, `ne_zero_of_secularDiag_eq_zero` and `pos_of_secularDiag_eq_zero`
(`Secular.lean:302`, `:308`), prove this. Each of the four theorems now derives `hq` and `hnn`
from the root hypothesis `e` with `have`. The block forms `lamMax_Rstack_eq_gammaTop` and
`topSpace_Rstack_eq_span` no longer build `hq` and `hnn` either: they pass the root hypothesis
straight through.

**Modeling choices** (`notes/archive/agent_reports/polish_frobenius_secular.md`, D21 of `notes/FLAGGED.md`).

1. The lemma is proved as finite matrix algebra (class B), not as a limit statement. The
   Bunch-Nielsen-Sorensen determinant identity is the route.
2. The scalar root `IsGammaTop` and the spectral statement `lamMax_Rmat_eq` are separate, so
   that the scalar optimization of `thm:stacksvd_weighted` does not depend on the matrix
   layer.

**Proof sketch.** The file has a scalar half and a spectral half, and they meet at
`secularDiag_stack`. The scalar half is in namespace `Scalars` of `StackSVDWeighted.lean`:
`secular` (`StackSVDWeighted.lean:97`) is the paper's `f`, `IsGammaTop`
(`StackSVDWeighted.lean:107`) is "root above `max_i w_i²`", `secular_strictMonoOn`
(`StackSVDWeighted.lean:136`) proves the strict monotonicity that the paper obtains by citing
the interlacing theorem of `bunch1978rank`, and `existsUnique_gammaTop`
(`StackSVDWeighted.lean:179`) adds existence by the intermediate value theorem, from `f ≤ 0` just
above `max_i w_i²` and `f → 1` at `+∞`. `Scalars.sum_secular_eq_neg_one`
(`StackSVDWeighted.lean:274`) restates the root as `∑_j p_j = 1`, the form `thm:stacksvd_weighted`
consumes. In the spectral half every claim reduces to the finite-rank resolvent layer
`RMT/R4.lean`, which this project proves for a Hermitian `W₀` plus one rank-one `q qᵀ`; the
specialization to `W₀ = Matrix.diagonal σ` is `resolv_diagonal` (`Secular.lean:113`) and
`R4secular_diagonal` (`Secular.lean:125`). `det_Rmat_sub` (`Secular.lean:137`) is the paper's
determinant identity, from `R4.det_sub_smul_one` (`RMT/R4.lean:237`) and the rank-one column
update `R4.det_one_add_col` (`RMT/R4.lean:501`), and claim 2 `mem_spectrum_iff_secularDiag`
(`Secular.lean:153`) then divides out the nonzero product `∏_i (σ_i - λ)`. Claim 1 is geometric
rather than a count of roots: `mulVec_eq_of_const_on` (`Secular.lean:178`) shows that a vector
supported on a level set of `σ` and orthogonal to `q` is an eigenvector at that value, and
`card_sub_one_le_finrank_eigenspace` (`Secular.lean:197`) counts dimensions, `|S|` for the span
of the coordinate vectors, minus at most one for the single condition `⟪ q, x ⟫ = 0`. Claim 3 is
four instantiations of R4 above the top of the diagonal, which `lamMax_diagonal`
(`Secular.lean:267`) identifies with `max_i σ_i`: `root_above_unique` (`Secular.lean:317`) from
`R4.eigenvalue_above_unique` (`RMT/R4.lean:550`), `lamMax_Rmat_eq` (`Secular.lean:329`) from
`R4.lamMax_eq` (`RMT/R4.lean:630`), `topSimple_Rmat` (`Secular.lean:340`) from `R4.topSimple`
(`RMT/R4.lean:672`), and `topSpace_Rmat_eq_span` (`Secular.lean:351`) from `R4.topSpace_eq_span`
(`RMT/R4.lean:656`), where `resolv_diagonal` turns the resolvent applied to `q` into the entries
`q_i/(σ_i - λ)`, so the paper's `ξ₁ ∝ (Σ - γ₁ I)⁻¹ ũ₀` comes out literally. Four bridge lemmas move
this to the paper's block matrix: `secularDiag_stack` (`Secular.lean:399`) collapses the stacked
sum to `Scalars.secular`, `iSup_sigmaStack` (`Secular.lean:420`) matches `max_p σ_p` with
`Scalars.wSqMax`, `isGammaTop_iff_stack` (`Secular.lean:436`) makes the scalar and the spectral
root conditions one statement, and `u0Stack_ne_zero` (`Secular.lean:443`) shows `ũ₀ ≠ 0`
directly (since F30 the four claim-3 theorems derive it from the root instead);
`lamMax_Rstack_eq_gammaTop` (`Secular.lean:468`) and `topSpace_Rstack_eq_span`
(`Secular.lean:481`) are the instances. Mathlib supplies `spectrum_diagonal` and
`finrank_span_eq_card`, and nothing else enters, because the section has no measure and no
limit; there is no Gaussian form and no Layer 2, since `heteroLaw_of_gaussian` proves
`HeteroLaw.align` directly and does not read this file (D21).

**Scope.** Same as the paper. `Secular.lean` is not on the proof path of
`heteroLaw_of_gaussian`, which proves `HeteroLaw.align` directly (D21).

**Audit trail.** `notes/archive/audit_global_2026-08-31.md` finding 4;
`notes/archive/agent_reports/polish_frobenius_secular.md`.

## 5. Comparisons, estimation, and the appendix items

### 5.1 `thm:stacksvd_binary_optimal_svd_stack`

**Paper.**

```latex
\begin{prop} \label{thm:stacksvd_binary_optimal_svd_stack}
    \stacksvd with weights $w_i = \mathds{1}\{\beta_i > 0\}$ dominates optimally-weighted \svdstack when $c_i \le 1$ for all $i$, yielding strict improvement when $\beta_{2}>0$.
\end{prop}
```

**Lean** (`StackSVD/Weighted.lean`, since F28; the scalar half is `Scalars.svdstackOpt_le_binary` and the
strict `Scalars.svdstackOpt_lt_binary`):

```lean
theorem thm_stacksvd_binary_optimal_svd_stack_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hc1 : ∀ i, c i < (m.tbl i).θ ^ 4 → c i ≤ 1) (hdet : ∃ i, c i < (m.tbl i).θ ^ 4)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW
        (fun i => if i ∈ Finset.univ.filter (fun i => c i < (m.tbl i).θ ^ 4) then 1 else 0)
        N ω)
      (Scalars.binaryStackSVDLimit (Finset.univ.filter fun i => c i < (m.tbl i).θ ^ 4)
        (fun i => (m.tbl i).θ) c)
    ∧ svdstackLimitOpt (fun i => beta ((m.tbl i).θ) (c i))
        ≤ Scalars.binaryStackSVDLimit (Finset.univ.filter fun i => c i < (m.tbl i).θ ^ 4)
            (fun i => (m.tbl i).θ) c
```

The strict version is `thm_stacksvd_binary_optimal_svd_stack_gaussian_strict`
(`StrictFacades.lean`), which replaces `hdet` by the cardinality hypothesis
`hcard : 2 ≤ (Finset.univ.filter fun i => c i < (m.tbl i).θ ^ 4).card` (at least two
supercritical tables).

**Hypotheses in words.**

- `[∀ N, IsProbabilityMeasure (μ N)]`, `m`, `c`: measures, tables, aspect ratios.
- `hc : ∀ i, 0 < c i`: `eq:RMT_limit`.
- `hc1 : ∀ i, c i < (m.tbl i).θ ^ 4 → c i ≤ 1`: the paper's `c_i ≤ 1`, weakened to the kept
  (detectable) tables only, which is all the scalar inequality needs.
- `hdet : ∃ i, c i < (m.tbl i).θ ^ 4`: the detectable set is nonempty. The paper assumes this
  through `cor.2`.
- `hreg`, `hG`: the regime of every table and the joint Gaussian law.

**Hypotheses beyond the paper.** `[NeZero M]` is technical: the shared `v` is read off table
`0` and the kept tables are stacked, so `Fin M` must be nonempty for the statement to
elaborate. `[∀ N, IsProbabilityMeasure (μ N)]` is technical. `hc : ∀ i, 0 < c i` is
paper-implicit, from `eq:RMT_limit` (`main_paper.tex:262`). `hreg : ∀ i, (m.tbl i).Regime (c i)`
is paper-implicit, from `assum:main` (`main_paper.tex:271`). `hdet : ∃ i, c i < (m.tbl i).θ ^ 4`
is paper-implicit: the paper's weights `w_i = 1{β_i > 0}` are the indicator of that set, and
`cor.2` (`main_paper.tex:429`) reads a limit on a nonempty subset only; with `hdet` false the
weight vector is zero and the estimator is undefined. `hc1 : ∀ i, c i < (m.tbl i).θ ^ 4 → c i ≤ 1`
is not beyond the paper, since it is the paper's `c_i ≤ 1` restricted to the kept tables. `hG :
m.JointGaussianNoise` has two halves. The Gaussian half is restrictive against the four-moment
class of `assum:general_noise` (`main_paper.tex:246`), user decision Q3; what is lost is every
non-Gaussian noise law with four bounded moments. The independence half is paper-implicit,
since `assum:general_noise` makes the entries i.i.d. across the table index `i` as well. The
strict facade replaces `hdet` by `hcard : 2 ≤ (Finset.univ.filter fun i => c i < (m.tbl i).θ ^ 4).card`,
which is the paper's own `β₂ > 0`.

**Modeling choices** (`notes/archive/class_A.md`, `notes/archive/agent_reports/class_A_scalars.md`).

1. The comparison is a conjunction of a convergence statement and a scalar inequality, so the
   reader sees both what converges and what is compared.
2. The weights are the indicator of `{i : c_i < θ_i⁴}`, the strict detectable set of `cor.2`
   (see the scope note of section 3.7 and `notes/paper_edits.md` E6).

**Proof sketch.** The statement is a convergence claim and a scalar inequality, and the two
halves share no step. The convergence is `cor.2` read at the detectable set.
`MultiTableModel.restrict` (`StackSVDWeighted.lean:1303`) rebuilds the kept tables as a
smaller model, `heteroLaw_binary_of_singleTableLaw` (`StackSVDWeighted.lean:1370`) applies
`prop_stacksvd_general` to that sub-model and rewrites the limit as
`Scalars.binaryStackSVDLimit` through `Scalars.Lw_binary` (`StackSVDWeighted.lean:627`), and
`stackPerfW_binary_tendsto_gaussian` (`StackSVD/Weighted.lean:187`) is the `align` field of the
resulting law. The scalar inequality is `Scalars.svdstackOpt_le_binary` (`Scalars.lean:258`),
which follows the paper's own three steps. Step 1 rewrites each svdstack summand as
`(θ_i⁴ − c_i)/(θ_i² + c_i) = (θ_i² − c_i) − c_i(1 − c_i)/(θ_i² + c_i)`
(`Scalars.svdTerm_decomp`, `Scalars.lean:201`). Step 2 bounds the correction with the
Chebyshev-type inequality `∑ a_i ≤ (∑ y_j)(∑ a_i/y_i)` (`Scalars.sum_le_sum_mul_sum_div`,
`Scalars.lean:189`); the hypothesis `c_i ≤ 1` is what keeps every `a_i = c_i(1 − c_i)`
nonnegative. Step 3 combines the two into `S · (T + C) ≤ (T² − C) − ((∑ c)² − ∑ c²)`
(`Scalars.binary_core`, `Scalars.lean:230`) and closes with the shape lemma
`s/(s+1) ≤ A/(A+B)` (`Scalars.div_succ_le_div_add`, `Scalars.lean:57`). The strict half swaps
`∑ c_i² ≤ (∑ c_i)²` for its strict form at two or more kept tables
(`Scalars.sum_sq_lt_sq_sum`, `Scalars.lean:209`), giving `Scalars.svdstackOpt_lt_binary`
(`Scalars.lean:293`), which `thm_stacksvd_binary_optimal_svd_stack_gaussian_strict`
(`StrictFacades.lean:52`) pairs with the same convergence. Layer 2 discharges the single-table
law of the restricted stack in `heteroLaw_binary_of_gaussian` (`StackSVDWeighted.lean:1430`) by
one call to `SpikedModel.singleTableLaw_of_gaussian` (`RMT/Full.lean:49`). From Mathlib the
scalar layer uses only `Finset` sum algebra and `Finset.sum_sq_le_sq_sum_of_nonneg`.

**Scope.** `hc1` is weaker than the paper's hypothesis. Otherwise identical. The svdstack
side's convergence is `thm_svdstack_weighted_gaussian` of section 4.2.

**Audit trail.** `notes/archive/class_A.md`,
`notes/archive/agent_reports/class_A_scalars.md`, `notes/archive/audit_independent_D33_2026-09-02.md`
section 1.

### 5.2 `prop:dominance`

**Paper.**

```latex
\begin{thm}
    \label{prop:dominance}
    Optimally weighted \stacksvd dominates unweighted \stacksvd and optimally weighted \svdstack, providing strict improvement above its recovery threshold when $\theta_i^2/c_i$ is not constant across $i$, and when at least two $\theta_i$ are nonzero (respectively).
\end{thm}
```

**Paper, as amended (E8).** The proposition gains the same hypothesis as `thm:stacksvd_weighted`, since its non-strict clause fails at `θ ≡ 0`; the exact replacement is E8.2 in `notes/paper_edits.md`.

```latex
\begin{thm}
    \label{prop:dominance}
    Assume $\theta_i > 0$ for at least one $i$. Optimally weighted \stacksvd dominates unweighted \stacksvd and optimally weighted \svdstack, providing strict improvement above its recovery threshold when $\theta_i^2/c_i$ is not constant across $i$, and when at least two $\theta_i$ are nonzero (respectively).
\end{thm}
```

**Lean** (`RMT/Het/Sup.lean`; the scalar halves are `Scalars.svdstackOpt_le_stackSVDLimitW`
and `Scalars.stackSVDLimit_le_stackSVDLimitW`):

```lean
theorem prop_dominance_gaussian [NeZero M]
    (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hθ : ∃ i, (m.tbl i).θ ≠ 0) (hreg : ∀ i, (m.tbl i).Regime (c i))
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => m.stackPerfW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
      (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c)
    ∧ svdstackLimitOpt (fun i => beta ((m.tbl i).θ) (c i))
        ≤ Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c
    ∧ stackSVDLimit (fun i => (m.tbl i).θ) c
        ≤ Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) c
```

The strict clauses are `prop_dominance_gaussian_strict`, `_strict_svdstack`,
`_strict_unweighted` (`StrictFacades.lean`). All three carry `hthr : 1 < ∑ i, θ_i⁴/c_i`;
`_strict_svdstack` carries `htwo` (two nonzero `θ_i`) and `_strict_unweighted` carries `hnc`
(`θ_i²/c_i` not constant), in the paper's order of "(respectively)"; `_strict` carries both.
None carries `hθ`, which `Scalars.exists_ne_zero_of_thr` (`Scalars.lean:372`) derives from `hthr`.

**Hypotheses in words.**

- `m`, `c`: the tables and the aspect ratios.
- `hc : ∀ i, 0 < c i`: `eq:RMT_limit`.
- `hθ : ∃ i, (m.tbl i).θ ≠ 0`: at `θ = 0` the optimal weight vector is zero, the weighted
  stack is the zero matrix, and the claim is false. Not in the paper; in the paper as amended
  (E8); the same hypothesis as in section 4.1.
- `hreg`, `hG`: the regime of every table and the joint Gaussian law.
- The theorem takes no `[IsProbabilityMeasure]` instance: `JointGaussianNoise` already pins
  each `μ N`.

**Hypotheses beyond the paper.** `[NeZero M]` is technical: `stackPerfW` reads the shared `v`
off table `0`, so `Fin M` must be nonempty. `hc : ∀ i, 0 < c i` is paper-implicit, from
`eq:RMT_limit` (`main_paper.tex:262`), and `hreg : ∀ i, (m.tbl i).Regime (c i)` is
paper-implicit, from `assum:main` (`main_paper.tex:271`). `hθ : ∃ i, (m.tbl i).θ ≠ 0` is
restrictive (finding E8, as in 4.1): at `θ ≡ 0` the optimal weights are all zero, the weighted stack is the zero matrix,
and its top eigenvector is arbitrary, so the paper's estimator is not defined there. The three
strict facades do not carry `hθ`, because `Scalars.exists_ne_zero_of_thr` (`Scalars.lean:372`)
derives it from their threshold hypothesis. `hG : m.JointGaussianNoise` splits as in 5.1: the
Gaussian half is restrictive against the four-moment class of `assum:general_noise`
(`main_paper.tex:246`, user decision Q3), and the independence half is paper-implicit. The
statement takes no `[IsProbabilityMeasure]` binder. The strict facades add
`hthr : 1 < ∑_i θ_i⁴/c_i`, `htwo` and `hnc`, which are the paper's own three side conditions
(above the recovery threshold, at least two nonzero `θ_i`, and `θ_i²/c_i` not constant).

**Modeling choices** (`notes/archive/class_A.md`).

1. The two dominance claims are separate conjuncts, so a reader can cite one without the
   other. The strict clauses live in a separate file to keep the main statement short.

2. The paper's two side conditions sit on the strict clauses only, in the paper's order of
   "(respectively)": `_strict_svdstack` takes `htwo` (at least two `θ_i` nonzero) and
   `_strict_unweighted` takes `hnc` (`θ_i²/c_i` not constant); all three strict facades take
   `hthr : 1 < ∑ i, θ_i⁴/c_i`.

**Proof sketch.** The theorem is one convergence and two scalar inequalities, each proved
separately. The convergence is `thm_stacksvd_weighted_gaussian` (`RMT/Het/Sup.lean:292`) read
at `w = optWstack θ c`, that is section 4.1 at the optimal weights; no new probabilistic step
enters. Both inequalities pass through one comparison lemma. `Scalars.gW` (`Scalars.lean:334`)
is the paper's `f(x) + 1`, `Scalars.stackSVDLimitW` (`Scalars.lean:408`) is its unique root in
`(0,1)`, and `Scalars.le_stackSVDLimitW` (`Scalars.lean:442`) puts any `y` in `[0,1]` with
`gW θ c y ≥ 1` below that root, because `gW` is strictly decreasing
(`Scalars.gW_strictAntiOn`, `Scalars.lean:358`). The svdstack half
`Scalars.svdstackOpt_le_stackSVDLimitW` (`Scalars.lean:489`) applies it at `y = S/(S+1)`, and
gets `gW θ c y ≥ 1` termwise from the paper's `eq:dominance_termwise_ineq`, proved as
`Scalars.dominance_term` (`Scalars.lean:468`). The unweighted half
`Scalars.stackSVDLimit_le_stackSVDLimitW` (`Scalars.lean:536`) applies it at
`y = (T² − C)/(T² + T)` and gets `gW θ c y ≥ 1` from Cauchy-Schwarz in Engel form, Mathlib's
`Finset.sq_sum_div_le_sum_sq_div`, at the terms `(θ_i²)²/(c_i(T² + T) + θ_i²(T² − C))`. The
same Mathlib inequality supplies the threshold `T²/C ≤ ∑_i θ_i⁴/c_i` that the comparison
lemma needs. Below the threshold both claims are immediate, since `stackSVDLimitW` is zero by
definition and both left sides are zero. The strict clauses use
`Scalars.svdstackOpt_lt_stackSVDLimitW` (`Scalars.lean:958`) and
`Scalars.stackSVDLimit_lt_stackSVDLimitW` (`Scalars.lean:1007`) in place of the two
inequalities. Layer 2 discharges the `HeteroLaw` at the optimal weights inside
`thm_stacksvd_weighted_gaussian`, by `heteroLaw_of_gaussian` (`RMT/Het/Sup.lean:205`) with
`Scalars.optWstack_ne_zero` (`RMT/Het/Sup.lean:88`) for its nonzero-weight side condition.

**Scope.** Same as the paper, except the Gaussian noise (the `hG` paragraph above).

**Audit trail.** `notes/archive/class_A.md`, `notes/archive/audit_independent_D33_2026-09-02.md`
section 1.

### 5.3 `prop:binarystacksvd_inadmissable`

**Paper.**

```latex
\begin{prop}
    \label{prop:binarystacksvd_inadmissable}
    For any $\eps\in(0,1)$, there exists a problem instance of size $M= \lceil e^{-\gamma} \exp(2/\eps)\rceil$ where optimally weighted \stacksvd has asymptotic performance greater than $1-\eps$, while optimally binary-weighted \stacksvd and optimally weighted \svdstack both fall below their recovery thresholds.
    %This is $M= O(\exp(2/\eps))$ using classical Bachmann-Landau notation.
    This additionally implies that unweighted \stacksvd and \svdstack are also below their recovery thresholds.
\end{prop}
```

**Lean** (`RMT/Het/Sup.lean`; the scalar half is `Scalars.binarystacksvd_inadmissable` and
`_ceil` in `Scalars.lean`):

```lean
theorem stackPerfW_binary_inad_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (hθ : ∀ i, (m.tbl i).θ = 1) (S : Finset (Fin M))
    [NeZero S.card] (hreg : ∀ i, (m.tbl i).Regime (Scalars.inadC M i))
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω) 0

theorem stackPerfW_opt_inad_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (hθ : ∀ i, (m.tbl i).θ = 1)
    {ε : ℝ} (hε : 0 < ε) (hε1 : ε < 1)
    (hM : Real.exp (-Real.eulerMascheroniConstant) * Real.exp (2 / ε) ≤ M)
    (hreg : ∀ i, (m.tbl i).Regime (Scalars.inadC M i))
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ
        (fun N ω => m.stackPerfW
          (Scalars.optWstack (fun i => (m.tbl i).θ) (Scalars.inadC M)) N ω)
        (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) (Scalars.inadC M)) ∧
      1 - ε < Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) (Scalars.inadC M)
```

The paper's last sentence, the one about the **unweighted** methods, is covered by the scalar
theorem, which states all five claims at once:

```lean
theorem Scalars.binarystacksvd_inadmissable {ε : ℝ} (hε : 0 < ε) (hε1 : ε < 1) (M : ℕ)
    (hM : Real.exp (-Real.eulerMascheroniConstant) * Real.exp (2 / ε) ≤ M) :
    1 - ε < stackSVDLimitW (inadTheta M) (inadC M) ∧
      binaryStackSVDLimitMax (inadTheta M) (inadC M) = 0 ∧
      svdstackLimitOpt (fun i => beta (inadTheta M i) (inadC M i)) = 0 ∧
      stackSVDLimit (inadTheta M) (inadC M) = 0 ∧
      svdstackLimit (fun i => beta (inadTheta M i) (inadC M i)) = 0
```

`Scalars.binarystacksvd_inadmissable_ceil` states the same at the paper's own size
`M = ⌈e^{-γ} exp(2/ε)⌉`, as an existence claim over `M` (scalar limits only; the model is
supplied below).

**Hypotheses in words.**

- `[∀ N, IsProbabilityMeasure (μ N)]`, `m : MultiTableModel μ M n d`: measures and tables.
- `hθ : ∀ i, (m.tbl i).θ = 1`: the instance is explicit, `θ_i = 1` in every table.
- `S : Finset (Fin M)` with `[NeZero S.card]`: any nonempty subset of the tables.
- `hε : 0 < ε`, `hε1 : ε < 1`: the paper's `ε ∈ (0,1)`.
- `hM : exp(-γ) * exp(2/ε) ≤ M`: the paper's `M = ⌈e^{-γ} exp(2/ε)⌉`, as a lower bound.
- `hreg : ∀ i, (m.tbl i).Regime (Scalars.inadC M i)`: the regime at the harmonic family
  `c_i = 2i + 1`.
- `hG : m.JointGaussianNoise`: independent Gaussian tables.

**Hypotheses beyond the paper.** `[NeZero M]` and `[∀ N, IsProbabilityMeasure (μ N)]` are
technical. `hθ : ∀ i, (m.tbl i).θ = 1` and `hreg : ∀ i, (m.tbl i).Regime (Scalars.inadC M i)`
are the paper's own instance (`θ_i = θ_0 = 1`, `c_i = 2i − 1` at `main_paper.tex:666`), not
extra hypotheses; the Lean family reads `c_i = 2i + 1` only because `Fin M` counts from `0`.
`S : Finset (Fin M)` is a quantifier rather than a hypothesis, and `[NeZero S.card]` is
technical: on the empty subset the weight vector is zero and the estimator is undefined, which
is why the bundled theorem writes `S.Nonempty` instead. `hε : 0 < ε` and `hε1 : ε < 1` are the
paper's `ε ∈ (0,1)`. `hM : exp(-γ) * exp(2/ε) ≤ M` is weaker than the paper's
`M = ⌈e^{-γ} exp(2/ε)⌉`, and the existential theorem states the ceiling itself. `hG :
m.JointGaussianNoise` splits as in 5.1: Gaussian noise is restrictive against the four-moment
class of `assum:general_noise` (`main_paper.tex:246`, user decision Q3), and independence
across tables is paper-implicit. The two general svdstack theorems quoted above add
`hc : ∀ i, 0 < c i` (paper-implicit, `eq:RMT_limit` at `main_paper.tex:262`), `hβdef` and
`hβ0` (derivable at the instance, in two lines, by `beta_one_inadC_eq_zero` of
`SVDStack/Inad.lean:27`) and `hwne : ∃ k, w k ≠ 0` (technical, and discharged at the instance
because `optW 0 = 1`). `prop_binarystacksvd_inadmissable_exists` keeps only `hε` and `hε1`.

**Modeling choices** (`notes/archive/class_A.md`, `notes/archive/agent_reports/class_A_rest.md`).

1. The paper's existence claim is stated twice: as a named instance (`Scalars.inadTheta`,
   `Scalars.inadC`) plus theorems about every model of that instance, and in the paper's own
   existential form with the model built (`Existence.lean`, below).
2. The binary statement quantifies over every nonempty `S`, so it covers the optimal binary
   weighting and the unweighted stack at once.
3. The harmonic bound uses `H_M > log M + γ` from Mathlib's Euler-Mascheroni constant.

**Proof sketch.** The proposition is five claims on one instance, and each one is an already
proved limit theorem read at `θ_i = 1`, `c_i = 2i + 1`. Three of the five are zero limits and
share one arithmetic fact: `(∑_{i ∈ S} θ_i²)² = |S|² ≤ ∑_{i ∈ S} c_i` for every subset,
because the first `|S|` odd numbers sum to `|S|²` and `S` selects `|S|` of them
(`Scalars.sq_card_le_sum_inadC`, `Scalars.lean:1316`). So the guard of
`Scalars.binaryStackSVDLimit` fails on every `S` (`Scalars.inad_binaryStackSVDLimit_eq_zero`,
`Scalars.lean:1331`), which settles optimally binary weighted stacksvd (the subset maximum)
and unweighted stacksvd (the subset `univ`) at once. The two svdstack claims follow from
`θ_i⁴ = 1 ≤ c_i`, so every `β_i` is zero (`beta_one_inadC_eq_zero`, `SVDStack/Inad.lean:27`)
and `Sval β = 0`. The fourth claim is the paper's harmonic estimate: `gW` of section 5.2 at
`x = 1 − ε` is at least `(ε/2) H_M` termwise, `H_M > log M + γ` comes from Mathlib's
`Real.eulerMascheroniConstant_lt_eulerMascheroniSeq'`, and `M ≥ e^{2/ε − γ}` then gives
`gW θ c (1 − ε) > 1` (`Scalars.one_lt_inad_gW`, `Scalars.lean:1421`), so the root exceeds
`1 − ε` by `Scalars.lt_stackSVDLimitW` (`Scalars.lean:914`).
`Scalars.binarystacksvd_inadmissable` (`Scalars.lean:1473`) collects the five scalar values and
`Scalars.binarystacksvd_inadmissable_ceil` (`Scalars.lean:1502`) restates them at the paper's
ceiling. On the model side each value is transported by its own limit theorem:
`stackPerfW_binary_inad_gaussian` (`RMT/Het/Sup.lean:440`) and `stackPerfW_opt_inad_gaussian`
(`RMT/Het/Sup.lean:456`) for the two stacksvd claims, `svdstackPerf_inad_gaussian`
(`SVDStack/Inad.lean:42`) and `svdstackPerfW_opt_inad_gaussian` (`SVDStack/Inad.lean:54`) for
the two svdstack claims, and `prop_stacksvd_general_gaussian` for the unweighted stack. Layer
2 discharges `HeteroLaw` at the optimal weights by `heteroLaw_of_gaussian`
(`RMT/Het/Sup.lean:205`) and the binary law by `heteroLaw_binary_of_gaussian`
(`StackSVDWeighted.lean:1430`). The existential form builds the model: `Sat.Rank1.inadModel`
(`Existence.lean:102`) is `Sat.Rank1.model` (`Sat.lean:82`) at `k_i = 2i + 1`, `l = 1`, so
`n_i N = (2i + 1)(N + 1)` and `d N = N + 1` realize `c_i = 2i + 1` exactly, and
`prop_binarystacksvd_inadmissable_exists` (`Existence.lean:124`) feeds it to
`prop_binarystacksvd_inadmissable_gaussian` (`Existence.lean:62`).

**Scope.** The instance is named, and the binary claim is uniform in `S`, so the theorems
about a model of the instance cover the paper's claim once a model exists; the existential
theorem below supplies one. The paper's closing sentence about the unweighted methods is the fourth and fifth
conjunct of `Scalars.binarystacksvd_inadmissable`; the probabilistic twin of the fourth is
`stackPerfW_binary_inad_gaussian` at `S = univ`.

The SVDstack clauses (optimally weighted and unweighted SVDstack below threshold) have
estimator theorems of their own in the tree (external statement audit finding 4). They
are the all-subcritical cases of
`thm:svd_stack_general` and `thm:svdstack_weighted`:

```lean
theorem thm_svd_stack_general_zero_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)) (hβ0 : ∀ i, β i = 0)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) 0

theorem thm_svdstack_weighted_zero_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (c β w : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)) (hβ0 : ∀ i, β i = 0) (hwne : ∃ k, w k ≠ 0)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerfW w N ω) 0
```

(`SVDStack/Main.lean`, `SVDStack/Weighted.lean`; both write `[NeZero M]`; see section 0
item 2.) At the instance every `β_i = beta 1 (2i + 1)` is `0`, because `θ_i⁴ = 1 ≤ c_i = 2i + 1`
and `beta` is `0` at and below the threshold (the third and fifth conjuncts of
`Scalars.binarystacksvd_inadmissable` are the scalar limits that follow), and the paper's optimal SVDstack weights `optW β i = 1/√(1 − β_i²)` are all `1`
there, so `thm_svdstack_weighted_zero_gaussian` applies at `w = optW β` with `hwne`
satisfied; the weight degeneracy of `AUDIT_DOC.md` 7.9 concerns the stackSVD weights at
`θ ≡ 0`, and `θ = 1` here. The composition at the instance is in the tree
(`SVDStack/Inad.lean`, item L3 of the external statement audit; the two facades there write
`[NeZero M]` in their own signatures):

```lean
theorem beta_one_inadC_eq_zero (M : ℕ) (i : Fin M) : beta 1 (Scalars.inadC M i) = 0

theorem svdstackPerf_inad_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (hθ : ∀ i, (m.tbl i).θ = 1)
    (hreg : ∀ i, (m.tbl i).Regime (Scalars.inadC M i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) 0

theorem svdstackPerfW_opt_inad_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) (hθ : ∀ i, (m.tbl i).θ = 1)
    (hreg : ∀ i, (m.tbl i).Regime (Scalars.inadC M i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ
      (fun N ω => m.svdstackPerfW (optW fun i => beta (m.tbl i).θ (Scalars.inadC M i)) N ω)
      0
```

These two theorems and the two stacksvd theorems above share the hypotheses `hθ`, `hreg`,
`hG`, so the four clauses of the proposition on the model hold on one instance under one
hypothesis set. "Below the recovery threshold" is read as "performance tends to `0` in
probability" in all four (`notes/archive/L3_inad_svdstack.md`, choice 1).

**The paper's existential form, with the model built** (`Existence.lean`; external
rank-one audit finding 1). The four theorems above and
`prop_stacksvd_general_gaussian` take a model of the instance as a hypothesis. The Gaussian
construction `Sat.Rank1.model` of `Sat.lean` (`n_i N = k_i (N+1)`, `d N = l (N+1)`, so
`c_i = k_i / l` exactly) supplies one at `k_i = 2i + 1`, `l = 1`: `Sat.Rank1.inadModel M`.
`prop_binarystacksvd_inadmissable_gaussian` bundles the five limits on any model of the
instance, and `prop_binarystacksvd_inadmissable_exists` is the paper's sentence:

```lean
theorem prop_binarystacksvd_inadmissable_gaussian [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)] (m : MultiTableModel μ M n d)
    (hθ : ∀ i, (m.tbl i).θ = 1) {ε : ℝ} (hε : 0 < ε) (hε1 : ε < 1)
    (hM : Real.exp (-Real.eulerMascheroniConstant) * Real.exp (2 / ε) ≤ M)
    (hreg : ∀ i, (m.tbl i).Regime (Scalars.inadC M i)) (hG : m.JointGaussianNoise) :
    (TendstoInProb μ
        (fun N ω => m.stackPerfW
          (Scalars.optWstack (fun i => (m.tbl i).θ) (Scalars.inadC M)) N ω)
        (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) (Scalars.inadC M)) ∧
      1 - ε < Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) (Scalars.inadC M)) ∧
    (∀ S : Finset (Fin M), S.Nonempty →
      TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω) 0) ∧
    TendstoInProb μ
      (fun N ω => m.svdstackPerfW (optW fun i => beta (m.tbl i).θ (Scalars.inadC M i)) N ω)
      0 ∧
    TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N)) 0 ∧
    TendstoInProb μ (fun N ω => m.svdstackPerf N ω) 0

theorem prop_binarystacksvd_inadmissable_exists {ε : ℝ} (hε : 0 < ε) (hε1 : ε < 1) :
    ∃ (M : ℕ) (_ : NeZero M) (Ω : ℕ → Type) (_ : ∀ N, MeasurableSpace (Ω N))
      (μ : ∀ N, Measure (Ω N)) (_ : ∀ N, IsProbabilityMeasure (μ N)) (n : Fin M → ℕ → ℕ)
      (d : ℕ → ℕ) (m : MultiTableModel μ M n d),
      M = ⌈Real.exp (-Real.eulerMascheroniConstant) * Real.exp (2 / ε)⌉₊ ∧
      (∀ i, (m.tbl i).θ = 1) ∧ (∀ i, (m.tbl i).Regime (Scalars.inadC M i)) ∧
      m.JointGaussianNoise ∧
      (TendstoInProb μ
          (fun N ω => m.stackPerfW
            (Scalars.optWstack (fun i => (m.tbl i).θ) (Scalars.inadC M)) N ω)
          (Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) (Scalars.inadC M)) ∧
        1 - ε < Scalars.stackSVDLimitW (fun i => (m.tbl i).θ) (Scalars.inadC M)) ∧
      (∀ S : Finset (Fin M), S.Nonempty →
        TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω) 0) ∧
      TendstoInProb μ
        (fun N ω => m.svdstackPerfW (optW fun i => beta (m.tbl i).θ (Scalars.inadC M i)) N ω)
        0 ∧
      TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N)) 0 ∧
      TendstoInProb μ (fun N ω => m.svdstackPerf N ω) 0
```

(`Existence.lean`; the first is in namespace `MultiTableModel`, the second in
`Sat.Rank1`.) The conjuncts, in order: the size is the paper's ceiling; the model has
`θ_i = 1`, the regime `c_i = 2i + 1`, and independent Gaussian noise; optimally weighted
stacksvd tends to a limit above `1 - ε`; every binary weighting of stacksvd tends to `0`;
optimally weighted svdstack tends to `0`; unweighted stacksvd tends to `0`; unweighted
svdstack tends to `0`. The witness `NeZero M` is part of the existential because the
estimators are defined for `M ≠ 0`; it holds since `M ≥ e^{2/ε - γ} > 1`. No hypothesis on a
model remains.

**Audit trail.** `notes/archive/class_A.md`,
`notes/archive/agent_reports/class_A_rest.md`, `notes/archive/L3_inad_svdstack.md`.

### 5.4 `thm:theta_est`

**Paper.**

```latex
\begin{prop} \label{thm:theta_est}
    Consider two tables $X_1,X_2$ following \Cref{assum:general_noise,assum:main}, where $\theta_1^4 > c_1$.
    Then, $\hat{\theta}_2$ in \Cref{eq:theta_estimation} is a consistent estimator of $\theta_2$, i.e. $\hat{\theta}_2 \pto \theta_2$.
\end{prop}
```

**Lean** (`ThetaEst.lean`):

```lean
theorem thm_theta_est (m : MultiTableModel μ M n d) {i j : Fin M} {ci cj : ℝ}
    (hci : 0 ≤ ci) (hthr : ci < (m.tbl i).θ ^ 4) (law : (m.tbl i).SingleTableLaw ci)
    (est : m.ThetaEstLaw i j cj) :
    TendstoInProb μ (fun N ω => m.thetaHat i j ci cj N ω) (m.tbl j).θ

theorem thm_theta_est_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) {i j : Fin M} (hij : i ≠ j) {ci cj : ℝ} (hci : 0 < ci)
    (hthr : ci < (m.tbl i).θ ^ 4) (hG : m.JointGaussianNoise)
    (hregi : (m.tbl i).Regime ci) (hregj : (m.tbl j).Regime cj) :
    TendstoInProb μ (fun N ω => m.thetaHat i j ci cj N ω) (m.tbl j).θ
```

**Lean, general noise** (`ThetaEst.lean` for Stage 0 of `notes/NONGAUSSIAN_SCOPE.md`,
2026-09-08, and `General/Layer1.lean` for Stage 1, 2026-09-09; the definitions `NoiseLaw`,
`noiseMatrix` and `JointGeneralNoise` are quoted in section 2.6):

```lean
theorem thetaEstLaw_of_general (m : MultiTableModel μ M n d) {ν : Measure ℝ} (hν : NoiseLaw ν)
    {i j : Fin M} (hij : i ≠ j) (hG : m.JointGeneralNoise ν) {cj : ℝ}
    (hregj : (m.tbl j).Regime cj) : m.ThetaEstLaw i j cj

theorem thm_theta_est_general [∀ N, IsProbabilityMeasure (μ N)] (m : MultiTableModel μ M n d)
    {ν : Measure ℝ} (hν : NoiseLaw ν) {i j : Fin M} (hij : i ≠ j) {ci cj : ℝ} (hci : 0 ≤ ci)
    (hthr : ci < (m.tbl i).θ ^ 4) (law : (m.tbl i).SingleTableLaw ci)
    (hG : m.JointGeneralNoise ν) (hregj : (m.tbl j).Regime cj) :
    TendstoInProb μ (fun N ω => m.thetaHat i j ci cj N ω) (m.tbl j).θ

theorem thm_theta_est_of_general [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) {i j : Fin M} (hij : i ≠ j) {ci cj : ℝ} (hci : 0 < ci)
    (hthr : ci < (m.tbl i).θ ^ 4) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.JointGeneralNoise ν)
    (hregi : (m.tbl i).Regime ci) (hregj : (m.tbl j).Regime cj)
    (hedge : ∀ ε > 0, Tendsto (fun N => μ N
      {ω | lamMax ((m.tbl i).W0 N ω) ((m.tbl i).isHermitian_W0 N ω)
        ≤ bulkEdge ci + ε}) atTop (𝓝 1)) :
    TendstoInProb μ (fun N ω => m.thetaHat i j ci cj N ω) (m.tbl j).θ
```

The form without `hedge`, `thm_theta_est_of_moments` (`General/Layer1.lean`, 2026-09-10), has
the same binders minus the edge; it is `thm_theta_est_of_general` applied to
`lamMax_W0_edge_of_general` at table `i`:

```lean
theorem thm_theta_est_of_moments [∀ N, IsProbabilityMeasure (μ N)]
    (m : MultiTableModel μ M n d) {i j : Fin M} (hij : i ≠ j) {ci cj : ℝ} (hci : 0 < ci)
    (hthr : ci < (m.tbl i).θ ^ 4) {ν : Measure ℝ} (hν : NoiseLaw ν)
    (hac : ν ≪ (volume : Measure ℝ)) (hG : m.JointGeneralNoise ν)
    (hregi : (m.tbl i).Regime ci) (hregj : (m.tbl j).Regime cj) :
    TendstoInProb μ (fun N ω => m.thetaHat i j ci cj N ω) (m.tbl j).θ
```

**Hypotheses in words.** `m : MultiTableModel μ M n d` is the multi-table model; `i j : Fin M`
name the two tables, `ci` and `cj` their aspect ratios; `hthr` is the paper's `θ_1⁴ > c_1`. Layer
1: `law` is the single-table law of table `i` (section 2.1) and `est` the two cross-table
limits of section 2.6, the noise of table `j` against the top direction of table `i`.
Gaussian: `hij : i ≠ j` keeps the tables distinct, `hG : m.JointGaussianNoise` is the Gaussian
noise, `hregi` and `hregj` are the asymptotic regimes of the two tables. The two tables are
named by indices `i` and `j` of one `MultiTableModel`, so the shared `v` is automatic.
General noise: `hν : NoiseLaw ν` says that the fixed law `ν` on `ℝ` is a probability measure
with mean `0`, second moment `1` and a finite fourth moment; `hG : m.JointGeneralNoise ν` says
that the `M` noise matrices are independent with i.i.d. entries of law `ν`; `law` is again the
single-table law of table `i`, kept as a hypothesis.

`thm_theta_est_of_general` (`General/Layer1.lean`) discharges the one hypothesis that
`thm_theta_est_general` keeps, `law`, the single-table law of table `i`. It pays for that with
the two binders of section 3.1, `hac : ν ≪ (volume : Measure ℝ)` and `hedge`, and the edge is
needed for table `i` alone: table `j` enters only through `ThetaEstLaw`, which Stage 0 proves
at `NoiseLaw` with no edge and no density. Two binders come back for the same reason:
`hregi : (m.tbl i).Regime ci`, because the single-table law of table `i` is now proved rather
than assumed, and `hci`, which strengthens from `0 ≤ ci` to `0 < ci`. The proof is one line,
`m.thm_theta_est` fed by `SpikedModel.singleTableLaw_of_general` and
`m.thetaEstLaw_of_general`.

**Hypotheses beyond the paper.** `hij : i ≠ j` (Gaussian only) is paper-implicit, since the
paper's `X_1` and `X_2` are two different tables; it also removes the need for a `[NeZero M]`
binder. `hci : 0 ≤ ci` (Layer 1) and `hci : 0 < ci` (Gaussian) are paper-implicit, from
`eq:RMT_limit` (`main_paper.tex:262`); the Layer 1 form is the weaker of the two, because it
also allows `c_i = 0`. `law : (m.tbl i).SingleTableLaw ci` is paper-implicit in the Layer 1
form: it is `prop:single_table`, which the paper's proof applies to `X_1`.
`est : m.ThetaEstLaw i j cj` is paper-implicit in the same way: its two fields are the two
concentration steps that the paper's proof takes from `lem:noise_projection_concentration`
and the independence of the tables (section 2.6). `hregi` and `hregj` (Gaussian only) are
paper-implicit, from `assum:main` (`main_paper.tex:271`). `[∀ N, IsProbabilityMeasure (μ N)]`
(Gaussian only) is technical. `hG : m.JointGaussianNoise` (Gaussian only) is restrictive in
its Gaussian half, as in 3.4: the paper asks only for the four-moment class of
`assum:general_noise` (`main_paper.tex:246`, user decision Q3), and the Layer 1 form holds
for every noise law that satisfies the single-table law and the two cross-table limits. The
independence half of `hG` is paper-implicit. Until F31 (2026-09-08) the Layer 1 form took
`hG` itself, the one rank-one exception (D35 of `notes/FLAGGED.md`). The general form
`thm_theta_est_general` (2026-09-08) removes the Gaussian restriction from the two cross-table
limits: `hν` and `hG : m.JointGeneralNoise ν` together are `assum:general_noise` of the paper
(a fixed law for every `N`, mean `0`, variance `1`, finite fourth moment; the paper's "fourth
moment bounded" reads, for one fixed law, as "finite"). What stays restrictive there is
`law : (m.tbl i).SingleTableLaw ci`, the single-table law of table `i` as a hypothesis. Its
general-noise discharge is Stage 1 of `notes/NONGAUSSIAN_SCOPE.md`, landed 2026-09-09:
`thm_theta_est_of_general` above drops `law` and takes `hac`, `hregi` and the edge of table
`i` in its place (section 3.1 for those two binders beyond the paper). `hregi` is not
needed in the general form, because the single-table law of table `i` enters as a hypothesis
and not through its Gaussian discharge; `hci : 0 ≤ ci` is the Layer 1 form.

**Modeling choices** (`notes/archive/thm_theta_est.md`).

1. The estimator `thetaHat` is the paper's `eq:theta_estimation` written as a function of the
   data. The quadratic root `eq:theta_est_quadratic` is a lemma in the same file.
2. `thetaHat` reads the paper's `v̂_1` as `ThetaEst.topDir` (`ThetaEst.lean:240`), the top
   right singular direction of table `i`, and the conclusion uses `0 ≤ θ` of both tables,
   the model field `hθ` (`ThetaEst.lean:891, 893`).
3. Item P of the roadmap, `SpikedModel.noise_projection_tendsto`
   (`lem:noise_projection_concentration`), is proved in this file and not assumed.

**Proof sketch.** The proof is a continuous mapping argument on two observables, so no
perturbation bound enters. The estimator `thetaHat` (`ThetaEst.lean:819`) is the scalar
function `ThetaEst.thetaHatFn` (`ThetaEst.lean:111`) evaluated at the pair
`(σ₁²(X_i), ‖X_j v̂_i‖²)`. The first coordinate converges by the `lamMax` field of the
`SingleTableLaw` of table `i`. The second is `sqNormMulVec_X_tendsto` (`ThetaEst.lean:825`),
which expands `‖X_j v̂_i‖²` into three parts (`SpikedModel.sqNormMulVec_X`,
`ThetaEst.lean:526`): a signal term `θ_j² ⟨v, v̂_i⟩²` that converges to `θ_j² β_i²` by the
`align` field of the same law, a cross term that tends to zero (the `crossTerm` field of
`ThetaEstLaw`), and the noise projection `‖E_j v̂_i‖²` that tends to `c_j` (its `noiseProj`
field). The last two carry the probabilistic content; for Gaussian noise they are
`cross_term_tendsto` (`ThetaEst.lean:742`) and `noise_projection_topDir_tendsto`
(`ThetaEst.lean:683`). Both are Chebyshev bounds after a Fubini step, `measure_randomDir_le`
(`ThetaEst.lean:638`): the direction `v̂_i` is a measurable function of the noise of table `i`,
independent of table `j`, so a second-moment bound uniform over unit directions still applies.
The variance itself is `ThetaEst.integral_sq_normSqMulVec` (`ThetaEst.lean:410`) for a fixed
direction, which is item P, `SpikedModel.noise_projection_tendsto` (`ThetaEst.lean:568`),
proved here and not assumed. The scalar side inverts `rhoSq`: `ThetaEst.thetaHat1_rhoSq`
(`ThetaEst.lean:122`) shows that the larger root of `eq:theta_est_quadratic` returns `θ_i`
above the threshold, `ThetaEst.thetaHatFn_limit` (`ThetaEst.lean:185`) evaluates the estimator
at the limit point and returns `θ_j`, and `ThetaEst.continuousAt_thetaHatFn`
(`ThetaEst.lean:194`) supplies the continuity, where `hthr` keeps the branch test of `betaSq`
locally constant. Layer 2 discharges the single-table law of table `i` in
`thm_theta_est_gaussian` (`ThetaEst.lean:897`) by `SpikedModel.singleTableLaw_of_gaussian`
(`RMT/Full.lean:49`) and the cross-table limits by `thetaEstLaw_of_gaussian`
(`ThetaEst.lean:812`). Mathlib supplies the Gaussian integrals and Chebyshev; the rest of the
section is proved here.

The general-noise chain (`ThetaEst.lean:906` to `1119`) repeats the three Gaussian steps with
the law swapped. `measure_randomDir_le_general` (`ThetaEst.lean:919`) is the Fubini step for
an arbitrary family of table laws `L` (the Gaussian and the general case both instantiate it).
`noise_projection_topDir_tendsto_general` (`ThetaEst.lean:962`) is the Chebyshev bound with
the second moment `NoiseLaw.integral_sq_normSqMulVec_le` (`Prob/NoiseMoments.lean`), the
inequality `(ν₄ + 2) n_j / d²` in place of the Gaussian identity `varSq · n_j / d²`, where
`ν₄ = ∫ x⁴ dν`; the bound tends to `0` because `n_j / d² = (n_j / d) / d → c_j · 0`.
`cross_term_tendsto_general` (`ThetaEst.lean:1024`) uses the exact second moment
`NoiseLaw.integral_sq_dotProduct_mulVec` (`= 1` at unit vectors), so its bound `1 / (d ε²)`
is the Gaussian one. Both second moments come from the linear-form moments of
`Prob/LinFormMoments.lean` (`LinForm.integral_linForm_sq`,
`LinForm.integral_linForm_pow_four_le`, the fourth moment of `∑ a_l y_l` at most
`ν₄ ∑ a_l⁴ + 3 (∑ a_l²)²`) and the product-law lemmas of `Prob/NoiseLaw.lean`
(`PiLaw.Centered`, a copy of `R2.Centered` at a general product law). `thetaEstLaw_of_general`
(`ThetaEst.lean:1098`) bundles the two limits, an `example` (`ThetaEst.lean:1106`) checks that
it gives `thetaEstLaw_of_gaussian` back at `ν = gaussianReal 0 1` (`noiseLaw_gaussian`,
`Prob/NoiseLaw.lean:143`), and `thm_theta_est_general` (`ThetaEst.lean:1114`) is
`thm_theta_est` applied to it. No declaration above `ThetaEst.lean:906` changed.

**Scope.** `thm_theta_est_gaussian`: same as the paper, except the Gaussian noise (choice 1
of section 3.1). `thm_theta_est_general`: the paper's noise class `assum:general_noise` for
the two cross-table limits, with `prop:single_table` of table `i` kept as a hypothesis.
`thm_theta_est_of_general`: the paper's noise class throughout, with a Lebesgue density and
the upper edge of table `i` added, so it is the only one of the three that assumes no limit
law.

**Audit trail.** `notes/archive/thm_theta_est.md`; the general form,
`notes/archive/thm_theta_est_general.md`; the Stage 1 form,
`notes/archive/prop_single_table_general.md` section 2.

### 5.5 `app:wstacksvd_mle`, the MLE identity and its marginalization

**Paper** (`main_paper.tex:1630` to `:1700`; the appendix has no theorem environment, so the
statement below is a paraphrase): the marginal maximum likelihood estimate of `v`,
after marginalizing Gaussian `u_i`, is
`v̂_MLE = v_max(∑_i θ_i²/(c_i+θ_i²) X_iᵀ X_i)`, which is exactly weighted stackSVD.

**Paper, as amended (E5), no change to any environment.** E5.2 names the step at `main_paper.tex:1651` as a standard Gaussian marginalization, so the reader separates it from the computation that follows. The Lean `thm_wstacksvd_mle_marginal` proves it.

**Lean** (`MLE.lean`, `MLEConverse.lean`):

```lean
theorem mleLogLik_eq [NeZero M] (m : MultiTableModel μ M n d) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (N : ℕ) (ω : Ω N) (v : Fin (d N) → ℝ) (hv : v ⬝ᵥ v = 1) :
    m.mleLogLik c N ω v
      = m.mleConst c N ω
        + (d N : ℝ) / 2
          * (v ⬝ᵥ m.stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω *ᵥ v)

theorem mleLogLik_max_iff_mem_topSpace [NeZero M] (m : MultiTableModel μ M n d) (c : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (N : ℕ) (ω : Ω N) (x : EuclideanSpace ℝ (Fin (d N))) (hx : ‖x‖ = 1) :
    (∀ v : Fin (d N) → ℝ, v ⬝ᵥ v = 1 →
        m.mleLogLik c N ω v ≤ m.mleLogLik c N ω (WithLp.ofLp x))
      ↔ x ∈ topSpace (m.stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ) c) N ω)
          (m.isHermitian_stackGramW _ N ω)
```

**Lean, the marginalization** (`MLEMarginal/Defs.lean`, `MLEMarginal/Main.lean`; L5). The
random-effects model of the appendix is a separate probability object, `reJointLaw n d θ v`:
table `i` is `θ_i u_i vᵀ + Z_i / √d` with `u_i ~ N(0, I/n_i)`,
`Z_i ~ gaussianMatrix (n i) d`, and the tables are independent (`Measure.pi`). Its Lebesgue
density is `reDensity n d θ v X = ∏_i ∏_k gaussDensity d (Σ_i(v)) (X_i)_k` with
`Σ_i(v) = mleCov d (θ_i² / (n_i / d)) v`, and `reLogLik = log reDensity`.

```lean
theorem thm_wstacksvd_mle_marginal [NeZero M] (m : MultiTableModel μ M n d) (N : ℕ)
    (ω : Ω N) (x : EuclideanSpace ℝ (Fin (d N))) (hx : ‖x‖ = 1) :
    (∀ v : Fin (d N) → ℝ,
        reJointLaw (fun i => n i N) (d N) (fun i => (m.tbl i).θ) v
          = (Measure.pi fun i => lebesgueMatrix (n i N) (d N)).withDensity
              (fun X => ENNReal.ofReal
                (reDensity (fun i => n i N) (d N) (fun i => (m.tbl i).θ) v X)))
    ∧ (∀ v : Fin (d N) → ℝ,
        reLogLik (fun i => n i N) (d N) (fun i => (m.tbl i).θ) v (fun i => (m.tbl i).X N ω)
          = m.mleLogLik (fun i => (n i N : ℝ) / d N) N ω v
            - ((∑ i, (n i N : ℝ)) * d N / 2) * Real.log (2 * Real.pi))
    ∧ ((∀ v : Fin (d N) → ℝ, v ⬝ᵥ v = 1 →
          reLogLik (fun i => n i N) (d N) (fun i => (m.tbl i).θ) v
              (fun i => (m.tbl i).X N ω)
            ≤ reLogLik (fun i => n i N) (d N) (fun i => (m.tbl i).θ) (WithLp.ofLp x)
                (fun i => (m.tbl i).X N ω))
        ↔ x ∈ topSpace (m.stackGramW (Scalars.optWstack (fun i => (m.tbl i).θ)
              (fun i => (n i N : ℝ) / d N)) N ω) (m.isHermitian_stackGramW _ N ω))
```

Clause 1 is the marginalization: the law of the random-effects model has the density
`reDensity` (the rows of each table are independent `N(0, Σ_i(v))` vectors,
`reTableLaw_eq_withDensity`). Clause 2 says that `log reDensity` at the observed tables
is the paper's `ℓ(v)` (that is, `mleLogLik` at `c_i = n_i/d`) up to the constant
`-(∑_i n_i) d/2 · log(2π)`, which does not depend on `v`. Clause 3 is the converse theorem
`mleLogLik_max_iff_mem_topSpace` above, restated for `reLogLik`: the unit maximizers of the
marginal log-likelihood are exactly the top eigenspace of the weighted stack Gram matrix
with the weights `θ_i²/(c_i + θ_i²)`.

**Hypotheses in words.**

- `m : MultiTableModel μ M n d`, `c : Fin M → ℝ`: the tables and the aspect ratios.
- `hc : ∀ i, 0 < c i`: `c_i > 0`. It makes the covariance `Σ(θ_i²/c_i, v)` invertible.
- `N : ℕ` and `ω : Ω N`: the first two statements are deterministic, at a fixed size and a
  fixed sample point. There is no limit and no measure hypothesis.
- `hv : v ⬝ᵥ v = 1`, `hx : ‖x‖ = 1`: unit vectors, as the paper has.
- `thm_wstacksvd_mle_marginal` takes only `m`, `N`, `ω`, `x`, `hx`. The positivity of
  `n_i` and `d` that the row density needs (`0 < n`, `0 < d` in `reRowLaw_eq_withDensity`)
  comes from the model fields `(m.tbl i).hn N`, `(m.tbl 0).hd N`. The random-effects law is
  a separate measure; nothing is assumed about the law of `m`.

**Hypotheses beyond the paper.** `[NeZero M]` is technical: the sum runs over `Fin M` and the
marginalization reads `d N` off table `0`. `hc : ∀ i, 0 < c i` is paper-implicit, from
`eq:RMT_limit` (`main_paper.tex:262`); it is also what makes `Σ_i(v)` invertible, which the
appendix's Sherman-Morrison step needs. `hv : v ⬝ᵥ v = 1` and `hx : ‖x‖ = 1` are the paper's
own `‖v‖₂ = 1`. `N : ℕ` and `ω : Ω N` are the point at which the identity holds, not
hypotheses; there is no limit, no regime and no measure hypothesis, so `assum:main`
(`main_paper.tex:271`) never enters this section. `thm_wstacksvd_mle_marginal` takes nothing
further: the positivity `0 < n_i` and `0 < d` that the density needs comes from the model
fields `(m.tbl i).hn N` and `(m.tbl 0).hd N`. Gaussian noise is not a restriction here,
because `reJointLaw` is a separate measure that the appendix itself defines as Gaussian, and
the theorem says nothing about the law of `m`; `assum:general_noise` (`main_paper.tex:246`) is
therefore not in play. No hypothesis of the three statements beyond `[NeZero M]` narrows the
paper's claim.

**Modeling choices** (`notes/archive/mle_identity.md`; question Q15 item 1 and D34 of `notes/FLAGGED.md`).

1. `mleLogLik` is a **definition** in Lean, written from the paper's display, with `c_i` in
   place of `n_i/d`. The Gaussian marginalization that produces that display is also
   proved (`thm_wstacksvd_mle_marginal`, item L5 of the external statement
   audit, `notes/archive/L5_mle_marginal.md`): the random-effects model is the measure
   `reJointLaw`, its Lebesgue density is `reDensity`, and `log reDensity` at the data is
   `mleLogLik` at `c_i = n_i/d` up to the constant `(∑_i n_i) d/2 · log(2π)`. The split
   between the two halves remains visible in the Lean names, which is what `notes/paper_edits.md`
   E5 item 3 asks the paper to make visible as well.
2. The identity is stated as "constant plus `d/2` times the Rayleigh form", not as an
   `argmax`. The converse theorem turns it into the two-sided statement about maximizers.
3. `Matrix.det_one_add_replicateCol_mul_replicateRow` (the Weinstein-Aronszajn identity at a
   rank-one update) and Sherman-Morrison come from Mathlib. `MLE.lean` proves the
   Sherman-Morrison instance it needs as `inv_mleCov`.
4. The random-effects model uses `c_i = n_i/d` at finite `N`, as the paper's display does;
   the asymptotic `c_i` of `SpikedModel.Regime` does not enter. The rows of each table are
   independent `N(0, Σ_i(v))` vectors (`reTableLaw_eq_pi`, `reTableLaw_eq_withDensity`) and
   the tables are independent (`Measure.pi`). The covariance `Σ_i(v)` is positive definite
   for every `v`, unit or not (`det Σ_i(v) = d^{-d}(1 + (θ_i²/c_i)‖v‖²)`), so clauses 1
   and 2 hold for every `v`; the unit norm enters only in clause 3.
5. No Mathlib multivariate Gaussian is used (v4.33.0 has none with a density). The row law
   is identified with the image of the standard Gaussian under the explicit square root
   `sqrtMleCov d a v = (√d)⁻¹(I + b vvᵀ)`, `b = a/(1 + √(1 + a‖v‖²))`, by characteristic
   functions (`Measure.ext_of_charFun`), and the density of that image comes from the
   linear change of variables `Measure.map_linearMap_addHaar_pi_eq_smul_addHaar`
   (`Prob/GaussianDensity.lean`, `Prob/WithDensityPi.lean`).

**Proof sketch.** The identity is algebra at a fixed `N` and `ω`. Two rank-one facts carry
it. First, `det Σ_i(v) = d^{-d}(1 + θ_i²/c_i)` at a unit `v` (`det_mleCov`, `MLE.lean:60`) does
not see the direction of `v`, so it joins the constant. Second, Sherman-Morrison
`Σ_i(v)⁻¹ = d(I − (a/(1+a)) vvᵀ)` (`inv_mleCov`, `MLE.lean:71`) turns `tr(Σ_i(v)⁻¹ X_iᵀX_i)`
into `d · tr(X_iᵀX_i) − d(a/(1+a)) vᵀX_iᵀX_i v` (`trace_inv_mleCov_mul`, `MLE.lean:84`). The
coefficient `a_i/(1 + a_i)` at `a_i = θ_i²/c_i` is the square of the paper's optimal weight
(`Scalars.optWstack_sq`, `MLE.lean:100`), which is why the Rayleigh form is the weighted stack
Gram matrix; `mleLogLik_eq` (`MLE.lean:161`) is then one `Finset` rearrangement. Mathlib
supplies the Weinstein-Aronszajn determinant
`Matrix.det_one_add_replicateCol_mul_replicateRow`; the Sherman-Morrison instance is proved
here. The converse `mem_topSpace_of_mleLogLik_max` (`MLEConverse.lean:128`) turns a maximizer
into the equality case of the Rayleigh bound, and `inner_toOp_self_eq_lamMax_iff`
(`MLEConverse.lean:78`) reads that case off by expanding `x` in the sorted eigenbasis, where
equality forces `∑_i (λ_max − λ_i) ⟨e_i, x⟩² = 0` with every term nonnegative. The
marginalization is the longer half, in four steps. One row of table `i` is `θ_i u_i v + z/√d`,
and `reRowLaw_eq_map_sqrtMleCov` (`MLEMarginal/RowLaw.lean:160`) identifies its law with the
image of the standard Gaussian under the explicit square root `sqrtMleCov`
(`MLEMarginal/Defs.lean:46`), by characteristic functions, since Mathlib v4.33.0 has no
multivariate Gaussian with a density. The density of that image is a linear change of
variables (`Prob/WithDensityPi.lean:93`, wrapped as `map_mulVec_pi_gaussianReal` in
`Prob/GaussianDensity.lean:122`). The rows recombine into a table by `reTableLaw_eq_pi`
(`MLEMarginal/TableLaw.lean:49`) and `reTableLaw_eq_withDensity`
(`MLEMarginal/TableLaw.lean:133`), and the tables into the product measure by
`reJointLaw_eq_withDensity` (`MLEMarginal/TableLaw.lean:191`). Taking logs, `reLogLik_eq`
(`MLEMarginal/LogLik.lean:62`) matches `log reDensity` with `mleLogLik` at `c_i = n_i/d` up to
`(∑_i n_i) d/2 · log(2π)`, and `thm_wstacksvd_mle_marginal` (`MLEMarginal/Main.lean:58`)
bundles the three clauses. There is no Layer 2 step, since the section takes no random matrix
theory hypothesis.

**Scope.** The Lean statements are the finite-`N` identity, its converse, and the
marginalization from the random-effects model to the written objective. Nothing asymptotic
is claimed, and nothing about the distribution of `m`.

**Audit trail.** `notes/archive/mle_identity.md`,
`notes/archive/agent_reports/mle_identity.md`, `notes/archive/agent_reports/q2_mle_converse_remark_ic.md`,
`notes/archive/audit_independent_D33_2026-09-02.md` section 4; for the marginalization,
`notes/archive/L5_mle_marginal.md` (statements, modeling choices, the per-unit delegation
record, the numeric audit with seed `20260903`).

### 5.6 The two remarks

**Paper.**

```latex
\begin{remark}\label{remark:stack_outperform_svd}
Unweighted \stacksvd can outperform optimally weighted \svdstack, particularly when many tables exhibit weak signal strengths below the detection threshold.
\end{remark}

\begin{remark}
\label{remark:svd_outperform_stack}
    Unweighted \svdstack can outperform binary-weighted \stacksvd when tables have highly imbalanced signal strengths.
\end{remark}
```

**Paper, as amended (E6), no change to the environment.** E6.1 rewrites the paragraph on the improvement over binary-weighted stackSVD, inside the justification at `main_paper.tex:603`: `θ_3 = (c_3+1)^{1/4}` puts table 3 inside `S`, so binary stackSVD is the full stack, whose limit tends to 0 as `c_3` grows, against `6/7` for SVDstack (0.1413 and 0.8570 at `c_3 = 10⁴`). The paragraph's
"binary-weighted stackSVD" is the rule of `cor.2` (keep the tables above the threshold,
`main_paper.tex:612`), not the maximum over subsets.

**Lean** (`Remarks.lean`, `RemarksUniform.lean`). Each remark is an explicit instance of the
closed forms. Declarations: `remark_stack_outperform_svd_stack`,
`remark_stack_outperform_svd_svdstack`, `remark_stack_outperform_svd_svdstack_uniform`;
`remark_svd_outperform_stack_two`, `remark_svd_outperform_stack_pair_stack`,
`remark_svd_outperform_stack_pair_svdstack`,
`remark_svd_outperform_stack_two_svdstack_unweighted`,
`remark_svd_outperform_stack_three_stack`, `remark_svd_outperform_stack_three_svdstack`,
`remark_svd_outperform_stack_three_binary`, `remark_three_stack_tendsto_zero`.

Example:

```lean
theorem remark_svd_outperform_stack_two (m : MultiTableModel μ 2 n₂ d)
    (hθ : ∀ i, (m.tbl i).θ = ![Real.sqrt 5, 4] i)
    (hreg : ∀ i, (m.tbl i).Regime (![1, 38.4] i)) (hG : m.JointGaussianNoise) :
    (∀ i, beta ((m.tbl i).θ) (![1, 38.4] i) ^ 2 = 4 / 5)
    ∧ TendstoInProb μ (fun N ω => m.stackPerfW 1 N ω) (2008 / 2310)
    ∧ (∀ S : Finset (Fin 2), S.Nonempty →
        TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω)
          (Scalars.binaryStackSVDLimit S (fun i => (m.tbl i).θ) ![1, 38.4])
        ∧ Scalars.binaryStackSVDLimit S (fun i => (m.tbl i).θ) ![1, 38.4] ≤ 2008 / 2310)
    ∧ TendstoInProb μ (fun N ω =>
        m.svdstackPerfW (optW fun i => beta ((m.tbl i).θ) (![1, 38.4] i)) N ω) (8 / 9)
    ∧ (2008 : ℝ) / 2310 < 8 / 9
```

**The existential form, with the model built** (`Existence.lean`; external
rank-one audit finding 2). "Can outperform" is an existence claim, and the
theorems above take the model as a hypothesis. `Sat.Rank1.oneModel M` (`θ_i = c_i = 1`,
`n_i N = d N = N + 1`) and `Sat.Rank1.twoModel` (`θ = (√5, 4)`, `c = (1, 38.4)`, so
`n_1 N = 5(N+1)`, `n_2 N = 192(N+1)`, `d N = 5(N+1)`) are the witnesses:

```lean
theorem remark_stack_outperform_svd_exists (M : ℕ) [NeZero M] (hM : 2 ≤ M) :
    ∃ (Ω : ℕ → Type) (_ : ∀ N, MeasurableSpace (Ω N)) (μ : ∀ N, Measure (Ω N))
      (_ : ∀ N, IsProbabilityMeasure (μ N)) (n : Fin M → ℕ → ℕ) (d : ℕ → ℕ)
      (m : MultiTableModel μ M n d),
      (∀ i, (m.tbl i).θ = 1) ∧ (∀ i, (m.tbl i).Regime 1) ∧ m.JointGaussianNoise ∧
      TendstoInProb μ (fun N ω => overlap (m.stack.X N ω) (m.stack.v N))
        (1 - 2 / ((M : ℝ) + 1)) ∧
      0 < 1 - 2 / ((M : ℝ) + 1) ∧
      TendstoInProb μ (fun N ω => m.svdstackPerf N ω) 0 ∧
      (∀ ε > 0, Tendsto (fun N => μ N {ω | ∃ w : Fin M → ℝ, (∃ k, w k ≠ 0) ∧
        ε ≤ m.svdstackPerfW w N ω}) atTop (𝓝 0))

theorem remark_svd_outperform_stack_exists :
    ∃ (Ω : ℕ → Type) (_ : ∀ N, MeasurableSpace (Ω N)) (μ : ∀ N, Measure (Ω N))
      (_ : ∀ N, IsProbabilityMeasure (μ N)) (n : Fin 2 → ℕ → ℕ) (d : ℕ → ℕ)
      (m : MultiTableModel μ 2 n d),
      (∀ i, (m.tbl i).θ = ![Real.sqrt 5, 4] i) ∧ (∀ i, (m.tbl i).Regime (![1, 38.4] i)) ∧
      m.JointGaussianNoise ∧
      (∀ i, beta ((m.tbl i).θ) (![1, 38.4] i) ^ 2 = 4 / 5) ∧
      (∀ S : Finset (Fin 2), S.Nonempty →
        TendstoInProb μ (fun N ω => m.stackPerfW (fun i => if i ∈ S then 1 else 0) N ω)
          (Scalars.binaryStackSVDLimit S (fun i => (m.tbl i).θ) ![1, 38.4]) ∧
        Scalars.binaryStackSVDLimit S (fun i => (m.tbl i).θ) ![1, 38.4] ≤ 2008 / 2310) ∧
      TendstoInProb μ (fun N ω =>
        m.svdstackPerfW (optW fun i => beta ((m.tbl i).θ) (![1, 38.4] i)) N ω) (8 / 9) ∧
      TendstoInProb μ (fun N ω => m.svdstackPerf N ω) (8 / 9) ∧
      (2008 : ℝ) / 2310 < 8 / 9
```

(both in namespace `Sat.Rank1`, `Existence.lean`). The first says: for every `M ≥ 2` there is
a Gaussian model of `M` tables at `θ_i = c_i = 1` on which unweighted stacksvd tends to
`1 - 2/(M+1) > 0`, unweighted svdstack tends to `0`, and with probability tending to one no
nonzero weighting of svdstack reaches any `ε > 0` (so optimally weighted svdstack is covered).
The `[NeZero M]` argument follows from `2 ≤ M` and is separate only because the estimators
are defined for `M ≠ 0`. The second says: there is a Gaussian model of two tables on which
every binary weighting of stacksvd (the unweighted stack included, at `S = univ`) tends to a
value at most `2008/2310 = 0.8693`, while optimally weighted and unweighted svdstack tend to
`8/9 = 0.8889`. The paper's `38.4` is `192/5`, which the row counts realize exactly.

**Hypotheses beyond the paper.** `[NeZero M]` and `[∀ N, IsProbabilityMeasure (μ N)]` are
technical. `hθ` and `hreg` are the paper's own instances, written as hypotheses on a model
rather than as a construction: `θ_i = c_i = 1` for `remark:stack_outperform_svd`
(`main_paper.tex:571`), and `θ = (√5, 4)`, `c = (1, 38.4)` for `remark:svd_outperform_stack`
(`main_paper.tex:617`). `hM : 2 ≤ M` in `remark_stack_outperform_svd_stack` is the paper's
`M > 1`. In `remark_svd_outperform_stack_pair_stack`, `hM : 4θ₀⁴ ≤ M c₀` is weaker than the
paper's `M > 4θ₀⁴/c₀`, since the equality case also gives zero; `hc : 0 < c₀` is
paper-implicit, from `eq:RMT_limit` (`main_paper.tex:262`); and `hne : i₀ ≠ i₁` is
paper-implicit, since the paper names two distinct tables. `hG : m.JointGaussianNoise` splits
as in 5.1: Gaussian noise is restrictive against the four-moment class of
`assum:general_noise` (`main_paper.tex:246`, user decision Q3), and independence across tables
is paper-implicit. The two existential theorems keep only `M` with `hM : 2 ≤ M` (the first)
and nothing at all (the second), because they build the model.

**Modeling choices** (`notes/archive/remark_facades.md`; question Q15 item 1 of `notes/FLAGGED.md`).

1. Each instance states the comparison with explicit rationals, so the claim is checkable
   without a plot.
2. The binary comparison is stated per subset instead of through
   `binaryStackSVDLimitMax`, so the maximum is visible.
3. The existential theorems fix one witness each (`oneModel M`, `twoModel`); the paper's
   other instances (two supercritical tables among `M - 2` null ones, and the `M = 3`
   example) stay conditional theorems in `Remarks.lean`.

**Proof sketch.** Each remark is an instance, so every proof evaluates closed forms and never
revisits the random matrix theory. For `remark:stack_outperform_svd` at `θ_i = c_i = 1`,
`thm_simple_thm1_stacksvd_gaussian` (`StackSVD/Main.lean:157`) gives the stacksvd limit
`1 − (1 + θ₀²)/(Mθ₀⁴ + θ₀²) = 1 − 2/(M + 1)`, whose guard `Mθ₀⁴ > θ₀²` holds at `M ≥ 2`
(`remark_stack_outperform_svd_stack`, `Remarks.lean:285`). Every table sits at its own
threshold `θ_i⁴ = c_i`, so `β = 0` and `thm_svd_stack_general_zero_gaussian`
(`SVDStack/Main.lean:398`) gives zero for unweighted svdstack (`Remarks.lean:303`). The
paper's stronger claim about optimally weighted svdstack is
`remark_stack_outperform_svd_svdstack_uniform` (`RemarksUniform.lean:60`), an instance at
`β = 0` of the Rayleigh bound `svdstackPerfW_uniform_bound` (`SVDStack/Rayleigh.lean:542`),
which is uniform over every nonzero weight vector, data dependent ones included. For
`remark:svd_outperform_stack` the two-table instance is arithmetic: both tables have
`β² = 4/5`, so `Sval β = 8` and optimally weighted svdstack is `8/9`
(`thm_svdstack_weighted_gaussian_opt`, `SVDStack/Weighted.lean:1418`), while
`Scalars.binaryStackSVDLimit` is `2008/2310` on `univ` and `4/5` on each singleton, and
`Finset (Fin 2)` has three nonempty subsets, which `decide` enumerates
(`remark_svd_outperform_stack_two`, `Remarks.lean:502`). Unweighted svdstack reaches the same
`8/9`, because equal `β_i` make `A_β` equicorrelated (`Scalars.svdstackLimit_const`,
`Scalars.lean:789`; `remark_svd_outperform_stack_two_svdstack_unweighted`,
`RemarksUniform.lean:89`). The `M`-table variant of that remark needs the top eigenpair of
`A_β` at `β = β₀(e_{i₀} + e_{i₁})`, namely `1 + β₀²` with eigenvector `(e_{i₀} + e_{i₁})/√2`,
computed in `lamMax_Abeta_pair` (`Remarks.lean:104`) and `svdstackLimit_pair`
(`Remarks.lean:248`). Layer 2 enters only inside the theorems just named, each of which
discharges its own law by `SpikedModel.singleTableLaw_of_gaussian` (`RMT/Full.lean:49`) or
`heteroLaw_binary_of_gaussian` (`StackSVDWeighted.lean:1430`). The existential forms plug in
`Sat.Rank1.oneModel` (`Existence.lean:160`) and `Sat.Rank1.twoModel` (`Existence.lean:223`),
both built by `Sat.Rank1.model` (`Sat.lean:82`) with integer row counts that realize `c`
exactly (`38.4 = 192/5`).

**Scope, and one paper finding.** The `M = 3` instance of `remark:svd_outperform_stack` uses
`θ = (2, 2, c₃^{1/4})`, so table 3 sits exactly on the threshold `θ₃⁴ = c₃`. Under the strict
detectable set of `cor.2` that table is discarded, binary stackSVD converges to
`62/72 = 0.861`, and svdstack converges to `6/7 = 0.857`. The paper's claim is the other way
around at that instance. `remark_svd_outperform_stack_three_binary` states the true
comparison. The fix `θ₃ = (c₃+1)^{1/4}` restores the paper's claim; this is
`notes/paper_edits.md` E6.

**Audit trail.** `notes/archive/remark_facades.md`,
`notes/archive/agent_reports/remark_facades.md`, `scripts/numeric/check_remark_examples.py` (deterministic,
PASS), `notes/archive/audit_independent_D33_2026-09-02.md` section 5.

## 6. Rank `r`: unaligned subspaces (Section 7) and the exactly aligned case (Appendix E)

Where each result lives in the paper. Sections 6.2 to 6.5 below are **Section 7** of the main
text (`sec:general_rankr_unaligned`, `main_paper.tex:749`). Section 6.1 is **Appendix D**
(`sec:unaligned_rank_r`, `:1941`), the proof appendix of that section. Sections 6.6 and 6.7
are **Appendix E** (`sec:rank_r`, `:2290`), the exactly aligned scenario; they are corollaries
stated only there, not in the main text. Section 6.8 is stated and proved entirely inside
**Appendix D** (`sec:insufficiency_of_global_weights_stacksvd`, `:2087`); the main text only
points to it (`:912`). Section 6.9 is the reverse: its statement is in the main text
(`:915`, in the weighted stackSVD subsection that starts at `:906`), and its proof is in
**Appendix D** (`sec:appendix_insufficiency_single_weights_stacksvd`, `:2169`). Section 6.8
is formalized at Layer 1 and Gaussian (Track G); section 6.9 is formalized at Layer 1
(`prop_singleweight_suboptimality_of_law`, under the hypothesis `hlaw`) and, since 2026-09-07,
unconditionally on its Gaussian witness (`prop_singleweight_suboptimality_gaussian`, F18b).

The model of Section 7 is `assum:unaligned`:

```latex
\begin{assum}[Unaligned subspaces]\label{assum:unaligned}
    The matrices $\{X_i \}$ are generated according to the following model:
    \begin{equation*}
        X_i = U_i \Theta_i (V R_i)^{\top} + E_i
    \end{equation*}
    where $R_i \in \mathbb{O}(r,r_i)$ and
    \begin{align*}
    \operatorname{Rank}\left(\sum_i R_i R_i^\top \right) = r
    \end{align*}
    Additionally, the noise conditions of \Cref{assum:general_noise} and RMT scaling limits of \Cref{eq:RMT_limit} are satisfied. $\Theta_i \in \R^{r_i \times r_i}$ is diagonal with positive entries and $U_i \in \mathbb{O}(n_i, r_i)$.
\end{assum}
```

In Lean this is `UnalignedModelR μ M n d r rk` (general `r_i`) and `UnalignedModel μ M n d r`
(`r_i = 1`). `assum:rank_r` is the exactly aligned case: `rk = alignedRk M r` and `R_i = 1`.

### 6.1 `lem:general_rank_delocalization` (Appendix D)

**Paper** (`main_paper.tex:1948`).

```latex
\begin{lemma}\label{lem:general_rank_delocalization}
\begin{equation*}
    \begin{split}
        \tilde{V} V &\ip \diag(\beta) R_{\text{stack}} := B_R \\
        \tilde{V} \tilde{V}^{\top} &\ip A_{\beta, R}.
    \end{split}
\end{equation*}
\end{lemma}
```

**Lean** (`RankR/GeneralMain.lean` at general `r_i`; `RankR/Unweighted.lean` at `r_i = 1`, since F28):

```lean
theorem lem_general_rank_delocalization_general (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise)
    {i i' : Fin M} (hii : i ≠ i') (j : Fin (rk i)) (j' : Fin (rk i')) :
    TendstoInProb μ (fun N ω => ⟪m.vhatG i j N ω, m.vhatG i' j' N ω⟫_ℝ)
      (beta ((m.tbl i).θ j) (c i) * beta ((m.tbl i').θ j') (c i') *
        ((m.R i)ᵀ * m.R i') j j')
```

**Hypotheses in words.**

- `m : UnalignedModelR μ M n d r rk`: `assum:unaligned` with per-table spike counts `rk i`.
- `c : Fin M → ℝ`: the limiting aspect ratios.
- `law : ∀ i, (m.tbl i).TableLawR (c i)`: the rank-`rk` single-table law per table.
- `hG : m.IndepNoise`: independence across tables, the paper's own hypothesis. It is weaker
  than `JointGaussianNoise`, which fixes the law as well; see modeling choice 2.
- `{i i' : Fin M}` and `hii : i ≠ i'`: the two spikes are in different tables.
- `j : Fin (rk i)`, `j' : Fin (rk i')`: the two spike indices.

**Modeling choices** (`notes/archive/rank_r_general.md`).

1. The statement is the **off-diagonal** half only. Inside one table the inner product is
   exactly `0` for `j ≠ j'` and exactly `1` for `j = j'`, because the rows of `Ṽ` from one
   table are orthonormal. That half needs no limit and no hypothesis.
2. The noise predicate is `IndepNoise`, not `JointGaussianNoise`, so the Layer 1 statement
   holds for any independent tables.
3. The two displays of the paper are `gramR` (the `Ṽ Ṽᵀ → A_{β,R}` half) and `VtV_tendsto`
   (the `Ṽ V → B_R` half), both entrywise.

**Scope.** Same as the paper.

**Audit trail.** `notes/archive/rank_r_general.md`,
`notes/archive/audit_rank_r_general_2026-08-31.md`,
`notes/archive/audit_independent_trackB_2026-09-02.md` (verdict: on track, 0 blockers).

### 6.2 `prop:general_rank_unweighted_svdstack`

**Paper.**

```latex
\begin{prop}\label{prop:general_rank_unweighted_svdstack}
Assume that for each data matrix $X_i$, $\Theta_i$ has $r_i$ distinct singular values. Further assume that $\lambda_r(A_{\beta,R}) - \lambda_{r+1}(A_{\beta,R}) > 0$, where $\lambda_{\tilde{r} + 1}(A_{\beta,R}) := -\infty$. Then the output of \svdstack satisfies:
    \begin{align*}
    \| \hat{V}_{\svdstack}^{\top} V \|_F^2 &\ip \left\| \Lambda^{-1/2} Q^{\top} \diag(\beta) R_{\text{stack}} \right\|_F^2.
    \end{align*}
\end{prop}
```

**Lean** (`RankR/GeneralGaussian.lean`):

```lean
theorem prop_general_rank_unweighted_svdstack_general_gaussian
    [∀ N, IsProbabilityMeasure (μ N)] (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hrr : r ≤ rtot rk)
    (hgap : TopGap (ABlock β m.R) (isHermitian_ABlock β m.R) r)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRG N ω) (limitRG β m.R)
```

The `r_i = 1` twin is `prop_general_rank_unweighted_svdstack_gaussian` (`RankR/Unweighted.lean`, since F28).
The paper-literal `‖V̂ᵀ V‖_F²` forms are
`prop_general_rank_unweighted_svdstack_general_frobenius` (with a frame that exists with
probability tending to one) and `_frobenius_eig` (at the canonical frame, no frame
hypothesis), in `RankR/GeneralFrob.lean`.

**Hypotheses in words.**

- `[∀ N, IsProbabilityMeasure (μ N)]`, `m`, `c`, `β`: measures, the model, the aspect ratios,
  the per-spike overlaps.
- `hc : ∀ i, 0 < c i`: `eq:RMT_limit`.
- `hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)`: names `β`.
- `hreg : ∀ i, (m.tbl i).Regime (c i)`: the regime of every table.
- `hrr : r ≤ rtot rk`: the paper's `r ≤ r̃`, where `r̃ = ∑_i r_i`. Without it the top-`r`
  spectral objects of `A_{β,R}` are not defined.
- `hgap : TopGap (ABlock β m.R) _ r`: the paper's `λ_r(A_{β,R}) - λ_{r+1}(A_{β,R}) > 0`.
  `TopGap A hA r` asks every eigenvalue of index below `r` to be strictly above every
  eigenvalue of index `r` or more. At `r = r̃` it holds vacuously, which is exactly the
  paper's convention `λ_{r̃+1} := -∞`.
- `hG : m.JointGaussianNoise`: independent Gaussian tables.
- The paper's "`Θ_i` has `r_i` distinct singular values" is a model field
  (`SpikedModelR.hθanti`), not a hypothesis of the theorem.

**Modeling choices** (`notes/archive/rank_r_defs.md`, D13 of `notes/FLAGGED.md`).

1. Performance is `tr(Bᵀ specInvTop A r B)`, the trace form, which equals the paper's
   `‖Λ^{-1/2} Qᵀ diag(β) R_stack‖_F²` when the top-`r` frame exists.
2. The `_frobenius_eig` variant fixes the canonical frame `(topEigMat, topEigVal)`, so no
   frame hypothesis is needed.

**Scope.** Same as the paper, plus the model's ordering of spikes inside a table, and
without the paper's `Rank(∑ R_i R_iᵀ) = r`, which the proof does not use (section 1, the
note under `UnalignedModelR`).

**Audit trail.** `notes/archive/rank_r_defs.md`, `notes/archive/audit_rank_r_2026-08-30.md`,
`notes/archive/audit_rank_r_t7_2026-08-30.md`,
`notes/archive/audit_independent_trackB_2026-09-02.md`,
`notes/archive/audit_independent_trackC_2026-09-02.md` (verdict: pass with findings, 0 must-fix).

### 6.3 `prop:stacksvd_subspace`

**Paper.**

```latex
\begin{prop} \label{prop:stacksvd_subspace}
Under \Cref{assum:unaligned}
    \begin{equation*}
        \left\| \Vstacksvd^{\top} V \right\|_F^2 \ip \sum_{j=1}^r
        \frac{\lambda_j^2(C) - \|c\|_1}{\lambda_j^2(C) + \lambda_j(C)} \mathds{1}\left\{ \lambda_j^2(C) > \|c\|_1\right\}.
    \end{equation*}
\end{prop}
```

**Paper, as amended (E9), no change to the environment.** E9.3 adds one sentence at `main_paper.tex:842`: at a tie across index `r` the top-`r` eigenspace carries no `d × r` frame, so `V̂_stacksvd` denotes any orthonormal basis of it on the event that the gap at index `r` is positive, whose probability tends to 1 (R8).

**Lean** (`RankR/SubspaceGaussian.lean`):

```lean
theorem prop_stacksvd_subspace_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ M n d r) (hG : m.JointGaussianNoise) (cc : Fin M → ℝ)
    (hreg : ∀ i, (m.tbl i).Regime (cc i)) (hc : 0 < ∑ i, cc i) :
    TendstoInProb μ (fun N ω => m.perfStackR N ω)
      (limitStackR (fun i => (m.tbl i).θ) m.R cc)

theorem prop_stacksvd_subspace_general_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (hG : m.JointGaussianNoise) (cc : Fin M → ℝ)
    (hreg : ∀ i, (m.tbl i).Regime (cc i)) (hc : 0 < ∑ i, cc i) :
    TendstoInProb μ (fun N ω => m.perfStackRG N ω)
      (limitStackRG (fun i => (m.tbl i).θ) m.R cc)
```

Layer 1 forms: `prop_stacksvd_subspace` and `prop_stacksvd_subspace_general`, both in
`RankR/SubspaceMain.lean` since F28 (2026-09-08; were `RankR/Subspace.lean` and
`RankR/SubspaceG.lean`). Paper-literal Frobenius forms at
`r_i = 1`: `prop_stacksvd_subspace_frobenius` and `_frobenius_eig` (`RankR/Frobenius.lean`),
discharged by `prop_stacksvd_subspace_frobenius_eig_gaussian` (`RankR/SubspaceGaussian.lean`).
At general `r_i` (F23, `RankR/SubspaceGaussian.lean` section 4b) the same three
forms are `prop_stacksvd_subspace_general_frobenius`, `_general_frobenius_eig` and

```lean
theorem prop_stacksvd_subspace_general_frobenius_eig_gaussian [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)] (m : UnalignedModelR μ M n d r rk)
    (hG : m.JointGaussianNoise) (cc : Fin M → ℝ)
    (hreg : ∀ i, (m.tbl i).Regime (cc i)) (hc : 0 < ∑ i, cc i) :
    TendstoInProb μ
      (fun N ω => frobSq ((topEigMat (m.isHermitian_stackGramG N ω) (m.r_le_d N))ᵀ * m.V N))
      (limitStackRG (fun i => (m.tbl i).θ) m.R cc)
```

whose conclusion is the paper's own `‖V̂_stacksvdᵀ V‖_F²` at the canonical top-`r` eigenframe
`topEigMat` of the stack Gram matrix (`frobSq` is the sum of the squared entries), with the
same hypotheses as `prop_stacksvd_subspace_general_gaussian`. The gap event that makes the
projector-sum surrogate `perfStackRG` equal to the Frobenius norm is proved almost surely at
every large `N` (`ae_topGap_stackGramG`, `topGap_stackGramG_whp_of_gaussian`).

**Hypotheses in words.** Only the model, the regime and the Gaussian noise. No eigengap of
`C`, no distinct-spike condition. That matches `main_paper.tex:842`.

- `[NeZero M]`: at least one table. Without it the stack has no rows.
- `[∀ N, IsProbabilityMeasure (μ N)]`: each `μ N` is a probability measure.
- `m : UnalignedModel μ M n d r` (or `UnalignedModelR μ M n d r rk` at general `r_i`):
  `assum:unaligned`.
- `hG : m.JointGaussianNoise`: independent Gaussian tables.
- `cc : Fin M → ℝ`: the per-table aspect ratios.
- `hreg : ∀ i, (m.tbl i).Regime (cc i)`: the regime of every table.
- `hc : 0 < ∑ i, cc i`: only the **total** `‖c‖₁` has to be positive, which is weaker than
  `∀ i, 0 < cc i`. It is what the one-matrix stack law needs.

**Modeling choices** (`notes/archive/prop_stacksvd_subspace.md`, D17 and D18 of `notes/FLAGGED.md`).

1. `SubspaceLaw` has one field, `align`, per spike of `C`, and no `topGap` field: a
   `topGap` field would be unsatisfiable when `∑_i n_i N < r ≤ d N`, and no theorem needs it
   (`notes/FLAGGED.md`, D18).
2. The limit `limitStackR θ R c` is the paper's sum, with `betaSq` supplying both branches.

**Scope, one narrowing and one widening.** The widening: the model has no field
`Rank(∑ R_i R_iᵀ) = r`, so the paper's rank condition is not assumed (section 1). The
narrowing: at general `r_i` each table is a `SpikedModelR`, whose fields
`hθnn` and `hθanti` force nonnegative and distinct `θ_ij` inside a table. Only `hθanti` narrows
the paper: `assum:unaligned` (`main_paper.tex:762`) says "`Θ_i ∈ R^{r_i × r_i}` is diagonal
with positive entries", and `hθnn` (`0 ≤ θ`, F8) is wider than that. The paper makes no
ordering assumption (`main_paper.tex:767`) and assumes distinct entries for svdstack only
(`:768`). So `prop_stacksvd_subspace_general_gaussian` is narrower than the paper at a table
with a repeated `θ_ij`. No stacksvd proof in this section reads either field, although nine
files elsewhere in the tree do read them (`notes/FLAGGED.md` item (4)), so an order-free
model would carry the same proof (Track A audit finding F1, `notes/FLAGGED.md` D-A1);
`notes/REVISION_LIST.md` R18 carries the per-field, per-result survival table.

**Audit trail.** `notes/archive/prop_stacksvd_subspace.md`,
`notes/archive/audit_stacksvd_subspace_2026-08-30.md`,
`notes/archive/audit_independent_trackA_2026-09-02.md` (verdict: pass with findings, 1 must-fix, a
scope documentation item, since addressed by the docstring at `SubspaceGaussian.lean`).

### 6.4 `thm:gen_rank_weight_svdstak`

**Paper.**

```latex
\begin{thm}\label{thm:gen_rank_weight_svdstak}
Suppose that $\beta_{ij} > 0$ for all $i,j$.
Then the optimal weights for \svdstack under \Cref{assum:unaligned} are $W\opt = \diag(1/\sqrt{1 - \beta_{ij}^2})$.
Specifically,
\begin{equation*}
    \left\| \Vsvdstack(W\opt)^{\top} V \right\|_F^2 \ip L\opt := r - \sum_{\ell = 1}^r \lambda_{\tilde{r} + 1 - \ell} (A_{\beta, R}^{-1/2} (W\opt)^{-2} A_{\beta,R}^{-1/2})
\end{equation*}
\end{thm}
```

**Paper, as amended (E2).** The hypothesis `β_ij > 0` is dropped, which is what the Lean theorem carries; the exact replacement is E2.2 in `notes/paper_edits.md`. E2.4 rewrites Step 2 of the proof (`main_paper.tex:2062`), which reads "`diag(β)` is invertible by assumption", so that the bound holds at every `N` and the attainment at `W⋆` holds at every rank of `B_R`.

```latex
\begin{thm}\label{thm:gen_rank_weight_svdstack}
% Then the optimal weights for \svdstack under \Cref{assum:unaligned} are $W\opt = \diag(1/\sqrt{1 - \beta_{ij}^2})$. 
Under \Cref{assum:unaligned}, \svdstack can be optimally weighted with $W\opt = \diag(1/\sqrt{1 - \beta_{ij}^2})$. 
Specifically,
\begin{equation*}
    \left\| \Vsvdstack(W\opt)^{\top} V \right\|_F^2 \ip L\opt := r - \sum_{\ell = 1}^r \lambda_{\tilde{r} + 1 - \ell} (A_{\beta, R}^{-1/2} (W\opt)^{-2} A_{\beta, R}^{-1/2})
\end{equation*}
\end{thm}
```

**Lean.** At `r_i = 1` (`RankR/WeightedMain.lean`, since F28; was `RankR/Weighted.lean`):

```lean
theorem thm_gen_rank_weight_svdstak (m : UnalignedModel μ M n d r) (c β : Fin M → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i, β i = beta (m.tbl i).θ (c i))
    (hrM : r ≤ M)
    (hrankB : (BR β m.R).rank = r)
    (law : ∀ i, (m.tbl i).SingleTableLaw (c i)) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRW (optWR β) N ω)
      (limitROpt β m.R (by simpa using hrM))
```

At general `r_i`, with Gaussian noise, the paper-literal corollary and the strongest form
(`RankR/GeneralGaussian.lean`):

```lean
theorem thm_gen_rank_weight_svdstak_general_r_paper_gaussian
    [∀ N, IsProbabilityMeasure (μ N)] (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ)
    (hc : ∀ i, 0 < c i) (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hr : 0 < r) (hβpos : ∀ i j, 0 < β i j)
    (hrank : (∑ i, m.R i * (m.R i)ᵀ).rank = r) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRGW (optWG β) N ω) (limitOptG β m.R _) ∧
      ∀ W, TopGap (ABlockW W β m.R) _ r → 0 < eigenvalue r-1 →
        TendstoInProb μ (fun N ω => m.perfRGW W N ω) (limitRGW W β m.R) ∧
          limitRGW W β m.R ≤ limitOptG β m.R _

theorem thm_gen_rank_weight_svdstak_general_r_full_gaussian
    [∀ N, IsProbabilityMeasure (μ N)] (m : UnalignedModelR μ M n d r rk)
    (c : Fin M → ℝ) (β : (i : Fin M) → Fin (rk i) → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hr : 0 < r) (hrr : r ≤ rtot rk)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRGW (optWG β) N ω) (limitOptG β m.R _) ∧
      (∀ W, TopGap … → 0 < … → TendstoInProb μ (fun N ω => m.perfRGW W N ω) (limitRGW W β m.R)
        ∧ limitRGW W β m.R ≤ limitOptG β m.R _) ∧
      ∀ ε > 0, Tendsto (fun N => μ N {ω | ∃ W, limitOptG β m.R _ + ε ≤ m.perfRGW W N ω})
        atTop (𝓝 0)
```

(The two blocks above are abridged at the elaborated proof terms inside the binders; the file
carries them in full.)

**Hypotheses in words.** The two Gaussian blocks above are abridged **only** inside the
conclusion, at the elaborated proof terms that `limitOptG` and the eigenvalue index take
(`limitOptG β m.R _` stands for `limitOptG β m.R (by simpa using hrr)` in `_full_gaussian`
and for `limitOptG β m.R (rank_le_card_of_rank_BBlock (rank_BBlock_of_paper m.R hβpos hrank))`
in `_paper_gaussian`; `TopGap …` and `0 < …` stand for
`TopGap (ABlockW W β m.R) (isHermitian_ABlockW W β m.R) r` and
`0 < (isHermitian_ABlockW W β m.R).eigenvalues₀ ⟨r - 1, _⟩`). **No binder is dropped**, and
every one is listed here. The performance `perfRGW W` is the trace form
`tr((WṼV)ᵀ specInvTop(WṼṼᵀWᵀ, r) (WṼV))` (`RankR/General.lean:537`); it equals
`‖V̂_svdstack(W)ᵀ V‖_F²` at a top-`r` eigenframe (`frobSq_vhatSvdstackGW`,
`RankR/GeneralFrob.lean:322`), and the theorem in that Frobenius form,
`thm_gen_rank_weight_svdstak_general_r_frobenius_eig` (`RankR/GeneralFrob.lean:440`), is
proved under `rank B_R = r` only.

Layer 1, `thm_gen_rank_weight_svdstak` (`r_i = 1`):

- `m : UnalignedModel μ M n d r`, `c β : Fin M → ℝ`: the model, the aspect ratios, the
  per-table overlaps.
- `hc : ∀ i, 0 < c i`: `eq:RMT_limit`.
- `hβdef : ∀ i, β i = beta (m.tbl i).θ (c i)`: names `β`.
- `hrM : r ≤ M`: at `r_i = 1` the paper's `r ≤ r̃` reads `r ≤ M`.
- `hrankB : (BR β m.R).rank = r`: replaces the paper's `β_ij > 0` plus
  `Rank(∑ R_i R_iᵀ) = r`. It is strictly weaker, and at `r = 1` it reads `∃ k, β_k ≠ 0`, the
  hypothesis of the rank-one `thm_svdstack_weighted`.
- `law : ∀ i, (m.tbl i).SingleTableLaw (c i)`, `hG : m.JointGaussianNoise`: the random matrix
  theory input, and cross-table independence.
- There is **no** `hr : 0 < r` here; see modeling choice 3.

Gaussian at general `r_i`, both blocks:

- `[∀ N, IsProbabilityMeasure (μ N)]`, `m : UnalignedModelR μ M n d r rk`,
  `c : Fin M → ℝ`, `β : (i : Fin M) → Fin (rk i) → ℝ`.
- `hc`, `hβdef`, `hreg : ∀ i, (m.tbl i).Regime (c i)`, `hG : m.JointGaussianNoise`: as above,
  per table.
- `hr : 0 < r`: the shared subspace is nontrivial, so the top-`r` eigenframe exists. Both
  blocks keep it.
- `_paper_gaussian` then takes `hβpos : ∀ i j, 0 < β i j` and
  `hrank : (∑ i, m.R i * (m.R i)ᵀ).rank = r`, which are the paper's own two hypotheses.
- `_full_gaussian` takes `hrr : r ≤ rtot rk` in their place, that is only the paper's
  `r ≤ r̃`, and is therefore the wider statement.
- The third conjunct of `_full_gaussian` is the upper bound uniform in `W`: with probability
  tending to one, no weight matrix beats `L⋆ + ε`. That removes the eigengap side condition
  from the paper's optimality claim; this is `notes/paper_edits.md` E2.

**Modeling choices** (`notes/archive/thm_gen_rank_weight_svdstak.md`, D15 and D16 of `notes/FLAGGED.md`).

1. `rank B_R = r` is the hypothesis of the eigenframe form; it is equivalent to the eigengap
   at `W⋆`. The weaker pair `∃ k, 0 < β k` plus `Rank(∑ R_i R_iᵀ) = r` does not give that gap:
   at `M = 3`, `r = 2`, `R = (e₁, e₂, e₁)`, `β = (0.5, 0, 0.5)` the gap at `W⋆` is `0`
   (D16, with the counterexample). The attainment itself needs no rank hypothesis:
   `thm_gen_rank_weight_svdstak_norank` (`RankR/WeightedUpper.lean:799`, `r_i = 1`) and
   `thm_gen_rank_weight_svdstak_general_r_norank` (`RankR/WeightedUpperG.lean:289`, general
   `r_i`) prove the convergence of the trace form at `W⋆` from `hc`, `hβdef`, `r ≤ r̃`, the law
   and the Gaussian or independent noise.
2. The optimality claim is a conjunct of `thm_gen_rank_weight_svdstak_max`
   (`RankR/WeightedMain.lean:173`) and of the quoted `_full_gaussian`, not of the first quoted
   signature `thm_gen_rank_weight_svdstak`, which is the convergence alone.
3. There is no `0 < r` hypothesis on the `r_i = 1` form: at `r = 0` both sides are `0`.

**Scope.** The Lean statement is wider than the paper's on the hypothesis, and the
uniform-in-`W` bound is stronger than the paper's conclusion. The paper's removal remark at
`main_paper.tex:901` (a component with `β_ij = 0` can be dropped) is **not** formalized and is
wrong in general: removal can break `Rank(∑ R_i R_iᵀ) = r`, after which `W⋆` has no eigengap
(question Q8, with a counterexample and a Monte Carlo rate).

**Audit trail.** `notes/archive/thm_gen_rank_weight_svdstak.md`,
`notes/archive/audit_rank_r_weighted_2026-08-30.md`,
`notes/archive/audit_rank_r_weighted_post_2026-08-30.md`,
`notes/archive/audit_independent_trackB_2026-09-02.md`.

### 6.5 The worked example of Section 7

**Paper**: `eq:psi_equation` (`main_paper.tex:848`) sets `M = 2`, `r = 2`, `r_i = 1`,
`R_1 = (1, 0)ᵀ`, `R_2 = (sin ψ, cos ψ)ᵀ`, and gives the two closed forms
`eq:perf_rankr_ex_svdstack` and `eq:perf_rankr_ex_stacksvd`.

**Lean** (`RankR/Example.lean`, `RankR/SubspaceGaussian.lean`):

```lean
theorem example_perfR_tendsto_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ 2 n d 2) (θ c ψ : ℝ) (hc : 0 < c)
    (hθ : ∀ i, (m.tbl i).θ = θ) (hR : m.R = Rex ψ)
    (hreg : ∀ i, (m.tbl i).Regime c) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfR N ω) (svdstackEx (beta θ c) (Real.sin ψ))

theorem example_perfStackR_tendsto_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModel μ 2 n d 2) (θ c ψ : ℝ) (hc : 0 < c)
    (hθ : ∀ i, (m.tbl i).θ = θ) (hR : m.R = Rex ψ)
    (hreg : ∀ i, (m.tbl i).Regime c) (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfStackR N ω) (stacksvdEx θ c (Real.sin ψ))
```

Both are unconditional under Gaussian noise. `RankR/Example.lean` also proves the kink
locations `eq:kink_location`.

**Hypotheses in words.** `m : UnalignedModel μ 2 n d 2` fixes `M = r = 2` with `r_i = 1`,
since `UnalignedModel` is the `r_i = 1` class (external statement
audit finding 16);
`hθ : ∀ i, (m.tbl i).θ = θ` gives both tables the same strength; `hR : m.R = Rex ψ` fixes the
paper's alignment `R_1 = (1,0)ᵀ`, `R_2 = (sin ψ, cos ψ)ᵀ`; `hc : 0 < c` and
`hreg : ∀ i, (m.tbl i).Regime c` are `eq:RMT_limit` at a shared `c`; `hG` is the joint
Gaussian law. The Layer 1 forms are `example_perfR_tendsto` (`RankR/Example.lean`, hypothesis
`law : ∀ i, (m.tbl i).SingleTableLaw c`) and `example_perfStackR_tendsto`
(`RankR/SubspaceMain.lean`, since F28; hypothesis `law : m.SubspaceLaw (2 * c)`).

**Scope.** The example is stated on `UnalignedModel`, the `r_i = 1` model, which has no
`hθanti` field at all (section 1.1). So the equal strengths `θ_1 = θ_2 = θ` of the paper's
figure are inside the class for every `θ`, and no ordering question arises here.

**Audit trail.** `notes/archive/audit_rank_r_t7_2026-08-30.md`,
`notes/archive/audit_independent_trackA_2026-09-02.md`.

### 6.6 `thm:rank_r_svdstack` (Appendix E)

**Paper** (`main_paper.tex:2401`).

```latex
\begin{cor} [Specialization of \Cref{thm:gen_rank_weight_svdstak}]\label{thm:rank_r_svdstack}
Under \Cref{assum:rank_r,assum:general_noise}, assuming that a) $\theta_{ij} \neq \theta_{ik}$ for all $i$, for all $j\neq k$, and b) $S_{j}\neq S_k$ for all $j\neq k$,
rank-$r$ weighted \svdstack satisfies:
\begin{align*}
    \left(v_j^\top \hat{v}_{j,\svdstack}\right)^2 &\pto \frac{S_j}{S_j+1} \quad \text{ for $j=1,2,\hdots,r$},\\
    \left\| V^\top \hat{V}_\svdstack \right\|_F^2 &\pto \sum_j \frac{S_j}{S_j+1}.
\end{align*}
\end{cor}
```

**Paper, as amended (E4).** The corollary does not call itself a specialization of `thm:gen_rank_weight_svdstak`, since that theorem labels no component; the exact replacement is E4.1 in `notes/paper_edits.md`.

```latex
\begin{cor}\label{thm:rank_r_svdstack} %[Generalization of \Cref{thm:rank_r_svdstack}]
Under \Cref{assum:rank_r,assum:general_noise}, assuming that a) $\theta_{ij} \neq \theta_{ik}$ for all $i$, for all $j\neq k$, and b) $S_{j}\neq S_k$ for all $j\neq k$,
rank-$r$ weighted \svdstack satisfies:
\begin{align*}
    \left(v_j^\top \hat{v}_{j,\svdstack}\right)^2 &\pto \frac{S_j}{S_j+1} \quad \text{ for $j=1,2,\hdots,r$},\\
    \left\| V^\top \hat{V}_\svdstack \right\|_F^2 &\pto \sum_j \frac{S_j}{S_j+1}.
\end{align*}
\end{cor}
```

**Lean.** Both clauses, unconditional under Gaussian noise
(`RankR/GeneralGaussian.lean`, `RankR/WeightedUpperG.lean`, `RankR/AlignedMain.lean`):

```lean
theorem thm_rank_r_svdstack_aggregate_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hG : m.JointGaussianNoise) :
    TendstoInProb μ (fun N ω => m.perfRGW (optWG β) N ω)
      (∑ j, Sagg β j / (Sagg β j + 1))

theorem thm_rank_r_svdstack_component_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hS : StrictAnti (Sagg β)) (hG : m.JointGaussianNoise) (j : Fin r) :
    TendstoInProb μ
      (fun N ω => (((m.vhatSvdstackGW (optWG β) N ω (topEigMat …) (topEigVal …))ᵀ * m.V N) j j) ^ 2)
      (Sagg β j / (Sagg β j + 1))
```

`thm_rank_r_svdstack_component_of_model_gaussian` drops `hS`: the proof needs only the
separation `SaggSep β` (`RankR/AlignedComponent.lean:162`), which `saggSep_of_model` (`:535`)
derives from the model's own ordering. The aggregate clause reads `perfRGW (optWG β)`, the
trace form; its Frobenius form at an eigenframe is proved under `rank B_R = r` only
(`thm_gen_rank_weight_svdstak_general_r_frobenius_eig`, `RankR/GeneralFrob.lean:440`).

**Hypotheses in words.** The component block is abridged only inside the conclusion: the two
`…` stand for the same pair of arguments in both places, `(m.isHermitian_gramWG (optWG β) N ω)`
and `(by simpa using le_rtot_alignedRk hM)`, so the frame reads
`topEigMat (m.isHermitian_gramWG (optWG β) N ω) (by simpa using le_rtot_alignedRk hM)` and
`topEigVal` the same. No binder is dropped.

- `[∀ N, IsProbabilityMeasure (μ N)]`: each `μ N` is a probability measure.
- `m : UnalignedModelR μ M n d r (alignedRk M r)`: the exactly aligned family of
  `assum:rank_r`, `r` spikes in every table.
- `c : Fin M → ℝ`, `β : Fin M → Fin r → ℝ`: the aspect ratios and the per-spike overlaps.
- `hc : ∀ i, 0 < c i`: `eq:RMT_limit`.
- `hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)`: names `β`.
- `hreg : ∀ i, (m.tbl i).Regime (c i)`: the regime of every table.
- `hR : ∀ i, m.R i = 1` is the paper's `R_i = I_r` of `assum:rank_r`.
- `hM : 0 < M` gives `r ≤ r̃` (through `le_rtot_alignedRk`, since `r̃ = M r`).
- `hG : m.JointGaussianNoise`: independent Gaussian tables.
- `hS : StrictAnti (Sagg β)` is the paper's hypothesis (b) `S_j ≠ S_k`, in the sorted form.
  Only the component clause takes it; `thm_rank_r_svdstack_component_of_model_gaussian` drops
  it, because the proof needs only `SaggSep β`, which the model's own ordering supplies
  (`saggSep_of_model`).
- `j : Fin r`: the component the clause is about.
- The paper's hypothesis (a) `θ_ij ≠ θ_ik` is a model field (`SpikedModelR.hθanti`).
- The aggregate clause needs no hypothesis on `S_j`: a component with `S_j = 0` contributes
  `0` to the sum.

**Modeling choices** (`notes/archive/rankr_D4_plan.md`, decision item 14 of `FLAGGED.md`).

1. `v̂_j` is column `j` of `V̂_svdstack(W⋆)` at the **canonical** top-`r` eigenframe
   `(topEigMat, topEigVal)`. `IsTopEigFrame` does not order its columns, so a component clause
   over an arbitrary frame is false. A general-frame variant adds the conjunct `lam k = λ_k`
   to the event of probability tending to one.
2. No `S_j > 0` hypothesis. The necessity scan (seed 20260904) shows the clause holds at
   `S_j = 0` with limit `0`: `(v_1ᵀ v̂_1)²` falls from 0.0277 at `d = 200` to 0.0011 at
   `d = 1600`. So `notes/paper_edits.md` E4 item 2, which asked for `S_j > 0`, is withdrawn
   (`notes/archive/paper_edits_amendments.md` amendment 1), and the two documents agree. A
   later change to either one has to move the other.

**Scope, and a paper finding.** The paper calls this corollary a specialization of
`thm:gen_rank_weight_svdstak`. It is not: the component clause is not a corollary of the
aggregate Frobenius statement, and it needed its own proof (`RankR/AlignedComponent.lean`,
`RankR/AlignedMain.lean`). This is `notes/paper_edits.md` E4. One narrowing, the same as in
section 6.7: with `hR : ∀ i, m.R i = 1` and the model field `hθanti : StrictAnti θ` in every
table, every table ranks the `r` global components in the same strict order (a common
component order). The paper's hypothesis (a) asks only that the `θ_ij` be distinct inside
each table, so an instance where two tables rank the components differently is outside the
class of both Gaussian theorems; `hS : StrictAnti (Sagg β)` is a relabeling, the per-table
`hθanti` is the restriction (external statement audit finding 7).

**Audit trail.** `notes/archive/rankr_D4_plan.md`,
`notes/archive/audit_independent_trackD_2026-09-02.md` (verdict: pass with findings, 0 must-fix, 9
minor).

**Beyond the paper (F2).** Asymptotic orthogonality of the svdstack columns
themselves, at the aligned model's optimal weights `W⋆`: the off-diagonal entries of
`V̂_svdstack(W⋆)ᵀ V` tend to `0` in probability. The paper states asymptotic orthogonality for
the stacksvd columns only (`main_paper.tex:926`, repeated at `:2318`).

**Lean** (`RankR/AlignedOrth.lean`):

```lean
theorem thm_rank_r_svdstack_offdiag_of_sep (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hsep : SaggSep β) (law : ∀ i, (m.tbl i).TableLawR (c i)) (hG : m.IndepNoise)
    {k j : Fin r} (hkj : k ≠ j) :
    TendstoInProb μ
      (fun N ω => (((m.vhatSvdstackGW (optWG β) N ω
          (topEigMat (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM))
          (topEigVal (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM)))ᵀ * m.V N) k j) ^ 2) 0

theorem thm_rank_r_svdstack_offdiag_of_model_gaussian [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (β : Fin M → Fin r → ℝ) (hc : ∀ i, 0 < c i)
    (hβdef : ∀ i j, β i j = beta ((m.tbl i).θ j) (c i))
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hR : ∀ i, m.R i = 1) (hM : 0 < M)
    (hG : m.JointGaussianNoise) {k j : Fin r} (hkj : k ≠ j) :
    TendstoInProb μ
      (fun N ω => (((m.vhatSvdstackGW (optWG β) N ω
          (topEigMat (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM))
          (topEigVal (m.isHermitian_gramWG (optWG β) N ω)
            (by simpa using le_rtot_alignedRk hM)))ᵀ * m.V N) k j) ^ 2) 0
```

The Frobenius corollary `thm_rank_r_svdstack_offdiag_frob_of_sep` sums the off-diagonal squares
over every pair of `Fin r × Fin r`; its Gaussian facade is
`thm_rank_r_svdstack_offdiag_frob_of_model_gaussian`.

**Hypotheses in words.** The Layer 1 form `thm_rank_r_svdstack_offdiag_of_sep` carries
`hsep : SaggSep β` (the separation of the sorted aggregates), `law : TableLawR` and
`hG : m.IndepNoise`; the Gaussian facade `_of_model_gaussian` has the binders of
`thm_rank_r_svdstack_component_of_model_gaussian` above, with `{k j : Fin r}` and
`hkj : k ≠ j` the two distinct components compared in place of the single `j : Fin r`, and
is the only separation-free form. No hypothesis `hS` and no hypothesis `0 < S_j`: the proof squeezes the off-diagonal entry against
the component clause at two distinct sorted eigenvector indices (Bessel), not against the
aggregate clause, so a component with `S_j = 0` is allowed, as in the component clause.

**Scope.** New result beyond the paper. It needs no top-`r` eigengap of `A⋆` (which would need
every `S_j > 0`, since that gap holds only when the limit matrix has one): the proof uses
the two-index Bessel squeeze, not the aggregate-clause route (`notes/archive/F_batch_2026-09-05.md`;
`notes/FOLLOWUP_LIST.md` item F2).

### 6.7 `thm:rank_r_stacksvd` (Appendix E)

**Paper** (`main_paper.tex:2337`).

```latex
\begin{cor}\label{thm:rank_r_stacksvd}
Under \Cref{assum:rank_r,assum:general_noise}, and assuming that $\tilde{\theta}_{jj} \neq \tilde{\theta}_{jk}$ for all weightings $j$, for all components $k\neq j$ (\Cref{eq:stacksvd_apptildTheta}),
rank-$r$ weighted \stacksvd satisfies:
\begin{align*}
    \left(v_j^\top \hat{v}_{j,\stacksvd}\right)^2 &\pto \gamma_j \quad \text{ for $j=1,2,\hdots,r$},\\
    \left\| V^\top \hat{V}_\stacksvd \right\|_F^2 &\pto \sum_j \gamma_j.
\end{align*}
\end{cor}
```

**Paper, as amended (E7.1, applied 2026-09-09).** The user's own edit of `main_paper.tex` fixes the corollary a different way than the proposed E7.1 replacement: instead of restricting to the ordered class, it states the hypothesis as `λ_jj ≠ λ_jk`, the secular roots E7.3's amended `alg:rank_r_stacksvd` ranks by, so the corollary holds for unordered `θ` too; see the "As applied" block of E7.1 in `notes/paper_edits.md`.

```latex
\begin{cor}\label{thm:rank_r_stacksvd}%[Generalization of \Cref{thm:rank_r_stacksvd}]
Under \Cref{assum:rank_r,assum:general_noise}, and assuming that $\lambda_{jj} \neq {\lambda}_{jk}$ for all weightings $j$, for all components $k\neq j$ (\Cref{eq:stacksvd_lambda}),
rank-$r$ weighted \stacksvd satisfies:
\begin{align*}
    \left(v_j^\top \hat{v}_{j,\stacksvd}\right)^2 &\pto \gamma_j \quad \text{ for $j=1,2,\hdots,r$},\\
    \left\| V^\top \hat{V}_\stacksvd \right\|_F^2 &\pto \sum_j \gamma_j.
\end{align*}
\end{cor}
```

**Paper, as amended (E3 and E7), no change to the environment.** E3.1 gives `γ_j` the value 0 below the threshold in `eq:stacksvd_gammak` (`main_paper.tex:2331`), which sits outside the corollary; E7.2 asks that the outliers `ρ_jk` of the detectable components be distinct, in the sentence at `:2368`; E7.3 replaces the rank of `θ̃_jj` by the rank of the outlier `ρ_jj` among those outliers (`ρ_jk` is the map `MPhet.rhoHet` of the largest secular root `λ_jk`, the quantity `ellSup` in `RankR/Het/Scalars.lean:307` ranks), in line `alg_line:stacksvd_rank` of `alg:rank_r_stacksvd`. The two edits change the two definitions the statement reads, `γ_j` below the threshold and `ℓ_j`, and leave the environment as it is.

**Lean.** Layer 1 (`RankR/StackMain.lean`) and Gaussian (`RankR/Het/Sup.lean`):

```lean
theorem thm_rank_r_stacksvd (m : UnalignedModelR μ M n d r (alignedRk M r))
    (c : Fin M → ℝ) (law : m.HeteroLawR c) :
    (∀ j : Fin r, TendstoInProb μ
        (fun N ω => ⟪m.vhatStackR c j N ω, m.colVecG N j⟫_ℝ ^ 2)
        (Scalars.gammaR m.thetaAligned c j)) ∧
      TendstoInProb μ (fun N ω => m.frobSqStackR c N ω)
        (∑ j : Fin r, Scalars.gammaR m.thetaAligned c j)

theorem thm_rank_r_stacksvd_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r (alignedRk M r)) (c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hR : ∀ i, m.R i = 1) (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (hnz : ∀ j : Fin r, ∃ i, 0 < (m.tbl i).θ j) :
    (∀ j : Fin r, TendstoInProb μ
        (fun N ω => ⟪m.vhatStackR c j N ω, m.colVecG N j⟫_ℝ ^ 2)
        (Scalars.gammaR m.thetaAligned c j)) ∧
      TendstoInProb μ (fun N ω => m.frobSqStackR c N ω)
        (∑ j : Fin r, Scalars.gammaR m.thetaAligned c j)
```

Split forms: `thm_rank_r_stacksvd_proj`, `_inner`, `_frobenius`, each with a `_gaussian` twin.

**Hypotheses in words.**

- Layer 1 takes only `law : m.HeteroLawR c`. `hc`, `hR` and the separation are not read by
  the implication; they are conditions for the law to hold (decision item 16 of `FLAGGED.md`).
- `hR : ∀ i, m.R i = 1` and the model class `alignedRk M r` are `assum:rank_r`.
- `hc` and `hreg` are `eq:RMT_limit`; `hG` is `assum:general_noise` in the Gaussian case.
- The paper's separation `θ̃_jj ≠ θ̃_jk` is **not** a hypothesis: `SpikedModelR.hθanti`
  orders the strengths strictly in every table, so it holds inside the class
  (`Scalars.rhoHet_lt_of_lt`) as soon as some table carries component `j`.
- `hnz : ∀ j, ∃ i, 0 < (m.tbl i).θ j` (F8) says that at least one table carries
  each component. When it fails at `j`, every weight `w_ij = θ_ij/√(θ_ij² + c_i)` is `0`, so
  `X_stack^(j) = 0` and the paper's own `θ̃_jj ≠ θ̃_jk` fails too (for `r ≥ 2`; at `r = 1` the
  side condition is empty, and the paper's `γ_1` is undefined, since its equation reads
  `0 = 1`); `hnz` is implied by the paper's hypothesis for `r ≥ 2` and is extra at `r = 1`,
  and the amended statement (E7.1) carries it. Tables with `θ_ij = 0` have weight
  `0`, and the proof drops them through
  the sub-model of `RankR/Het/Sub.lean` (`UnalignedModelR.sub`, `key_of_gaussian`).

**Modeling choices** (`notes/archive/rankr_D1_statement.md`, `notes/archive/rankr_TrackE_plan.md` v2,
decision items 15 and 17 of `FLAGGED.md`).

1. `γ_j := Scalars.gammaR`, total: the rank-one limit `stackSVDLimitW` read on column `j`,
   with value `0` at or below the threshold `∑_i θ_ij⁴/c_i ≤ 1`. The paper defines `γ_j` only
   above the threshold; this is `notes/paper_edits.md` E3.
2. `ℓ_j := Scalars.ellR`, 0-based and total. The paper's 1-based rank is `ellR + 1`.
3. The overlap is `overlapIdx` at index `ℓ_j`; the sign of `v̂_j` is fixed by
   `⟪v̂_j, v_j⟫ ≥ 0`.
4. The Frobenius clause is `∑_j ∑_k ⟪v̂_j, v_k⟫²`, which needs the field `crossProj`.
5. `simpleIdxJ` asks simplicity at the one index `ℓ_j` (`SimpleIdx`), not `SimpleSpec` at
   `ℓ_j + 1`. The stronger form ruled out ties that no theorem reads.
6. The core Gaussian theorem is stated at `Scalars.ellSup`, the rank among the supercritical
   outliers, with a second bulk theorem at every index at or above `numSup`. Inside the class
   `ellSup = ellR`.

**Scope, and paper finding E7.** The spectral index of the statement is
`Scalars.ellR m.thetaAligned c j`, the paper's `ℓ_j`, in every place it occurs: inside
`stackOverlapJ` and `vhatStackR` (`RankR/StackGamma.lean:346, :362`) and inside the
`crossProj` and `simpleIdxJ` fields of `HeteroLawR`. It equals the component index `j` only
through `ellR_thetaAligned` (`RankR/StackGamma.lean:323`), which takes two side conditions,
`hc : ∀ i, 0 < c i` and `hex : ∃ i, 0 < (m.tbl i).θ j`. The Layer 1 theorem `thm_rank_r_stacksvd` assumes
**neither**, so read its index as the paper's `ℓ_j`; the Gaussian form
`thm_rank_r_stacksvd_gaussian` has both, so there `ℓ_j = j`.

The model class carries `hθnn` and `hθanti` in every table plus `R_i = 1`, which is the
paper's own observation at `main_paper.tex:2369` ("if `θ_ij` follow the same ordering in each
table, then the ordering is preserved for each weighting, and so `ℓ_j = j`"). The paper's
worked example of `alg:rank_r_stacksvd` (`main_paper.tex:2359`), `θ_1 = (2, 1)` and
`θ_2 = (1, 10)`, is therefore **not** an `UnalignedModelR` with `R_i = 1`. The paper allows an
unordered `θ` and asks only `θ̃_jj ≠ θ̃_jk` (`:2306`, with the condition and the conclusion at
`:2368` and `:2369`). At unordered `θ` the paper's `ℓ_j` is not the outlier rank, and more
than generality is at stake: with the paper's own `ℓ_j` the first display is **false** on an
explicit `M = 2`, `r = 2` instance, where index 0 carries `0.571 ± 0.008` and the paper's index
1 carries `0.011 ± 0.003` at `d = 1200`. Inside the ordered class the display holds. The
proposed fix (outlier rank, or an ordering hypothesis) is `notes/paper_edits.md` E7. `Scalars.ellR`
is kept general so the statement does not change if the class is widened.

**Numeric check.** The wave-2b audit instantiated `thm_rank_r_stacksvd_gaussian` on a
concrete Gaussian model in Lean: `M = 2`, `r = 2`, `n_i N = d N = N + 3` (so `c = (1, 1)`),
`θ = [[3, 3/2], [5/2, 1]]`, `μ N = Measure.pi (fun i => gaussianMatrix (N+3) (N+3))`. Every
hypothesis discharges, the instance compiles, and `#print axioms` on it gives the three
standard axioms, so the hypothesis set is satisfiable and the theorem is not vacuous. Both
components are supercritical there (`∑ θ_i0⁴/c_i = 120.06`, `∑ θ_i1⁴/c_i = 6.06`), so both
`γ_j` are positive. `γ_j` matched `theory_pred.R` (`compute_x_star`) to `1.0e-13`:
`0.928594873932` and `0.614906970741`. Monte Carlo at `d = 500`, 20 draws, seed 20260903 gave
`0.92845 ± 0.00082` and `0.61013 ± 0.00729` with Frobenius `1.53912` against `∑ γ_j = 1.54350`.

**Audit trail.** `notes/archive/rankr_D1_statement.md`, `notes/archive/rankr_TrackE_plan.md`,
`notes/archive/audit_rankr_plan_E_2026-09-02.md` (accept with corrections, all applied),
`notes/archive/audit_independent_TrackE_wave1_2026-09-02.md` (pass, 0 must-fix),
`notes/archive/audit_independent_TrackE_wave2a_2026-09-02.md` (pass, 0 must-fix, 5 should-fix),
`notes/archive/audit_independent_TrackE_wave2b_2026-09-02.md` (pass, 0 must-fix, 8 should-fix and 1
scope note).

### 6.8 `prop:gen_rank_stacksvd_singleweight` (Appendix D)

**Paper** (`main_paper.tex:2104` to `:2121`).

```latex
\begin{assum}\label{assum:gen_rank_stacksvd_eig_sep}
    The top $r$ eigenvalues of $G$, denoted $\gamma_1, \ldots, \gamma_r$, are distinct and satisfy 
    \begin{equation}
        \sup_{\ell \in [r]} \sum_{i=1}^M \frac{c_i w_i^4}{(\gamma_{\ell} - w_i^2)^2} < 1
    \end{equation}
    and $\gamma_r > \sup_{i \in [M]} w_i^2$.
\end{assum}

\begin{prop}\label{prop:gen_rank_stacksvd_singleweight}
    Under Assumption \ref{assum:gen_rank_stacksvd_eig_sep}, for $\ell = 1,\hdots,r$ the matrix 
    \begin{equation}
        \sum_{i=1}^M \left( \frac{w_i^2}{\gamma_{\ell} - w_i^2} \right)R_i \Theta_i^2 R_i^{\top}.
    \end{equation}
    has a unit eigenvector $z_{\ell}$ with eigenvalue $1$. The performance of single-table weighted \stacksvd satisfies
    \begin{equation}
        \|\Vstacksvd(w)^{\top} V\|_F^2 \ip \sum_{\ell = 1}^r \frac{1 - \sum_{i=1}^M \frac{c_i w_i^4}{(\gamma_{\ell} - w_i^2)^2}}{\gamma_{\ell} z_{\ell}^{\top} \left[ \sum_{i=1}^M \left( \frac{w_i^2}{(w_i^2 - \gamma_{\ell})^2} \right)R_i \Theta_i^2 R_i^{\top}  \right] z_{\ell}}
    \end{equation}
\end{prop}
```

**Paper, as amended (E11), no change to the environment.** E11.1 adds one sentence after the proposition: the distinctness of the `γ_ℓ` carries the value of the formula, not only the proof, because at a tied root the `z_ℓ` are fixed only up to a rotation inside the tied plane and a sum of reciprocals is not invariant.

`G = A Aᵀ + Σ` is the covariance of the paper's `X_stack(w)` (`main_paper.tex:2087` to
`:2103`, the setup at the top of this section). `Σ = diag(w_1² I_{n_1}, …, w_M² I_{n_M})`, and
`A` stacks the rows `w_i U_i Θ_i R_iᵀ`. The model is `assum:unaligned`, given in full at the
top of section 6.

**Lean.** Layer 1 (`RankR/SingleWeight/`) and Gaussian (`RankR/SingleWeight/Het/`, Track G
of `notes/archive/trackG_plan.md`, 2026-09-05). The deterministic scalar layer
(`RankR/SingleWeight/Scalars.lean`):

```lean
/-- `S_i = R_i Θ_i² R_iᵀ`, the signal covariance of table `i` in the shared basis. It is the
matrix that carries `Θ_i` into the displays of `main_paper.tex:2115` and `:2119`. -/
noncomputable def sigMat {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (i : Fin M) :
    Matrix (Fin r) (Fin r) ℝ :=
  R i * Matrix.diagonal (fun j => θ i j ^ 2) * (R i)ᵀ

/-- `∑_i w_i²/(γ - w_i²) R_i Θ_i² R_iᵀ`, the matrix of `main_paper.tex:2115`. -/
noncomputable def secMat {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w : Fin M → ℝ) (γ : ℝ) :
    Matrix (Fin r) (Fin r) ℝ :=
  ∑ i, (w i ^ 2 / (γ - w i ^ 2)) • sigMat θ R i

/-- `∑_i w_i²/(γ - w_i²)² R_i Θ_i² R_iᵀ`, the matrix inside the denominator of
`main_paper.tex:2119`. It is `- d/dγ (secMat θ R w γ)`. The paper writes the denominator
`(w_i² - γ)²`, which is the same number. -/
noncomputable def secDerivMat {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w : Fin M → ℝ) (γ : ℝ) :
    Matrix (Fin r) (Fin r) ℝ :=
  ∑ i, (w i ^ 2 / (γ - w i ^ 2) ^ 2) • sigMat θ R i

/-- `γ` is a root of the rank-`r` secular equation above `max_i w_i²`
(`main_paper.tex:2115`). The rank-one twin is `Scalars.IsGammaTop`
(`StackSVDWeighted.lean:107`). -/
def IsSecularRoot {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w : Fin M → ℝ) (γ : ℝ) : Prop :=
  Scalars.wSqMax w < γ ∧ (1 - secMat θ R w γ).det = 0

/-- The first half of **`prop:gen_rank_stacksvd_singleweight`** (`main_paper.tex:2112`): a
root of the secular equation gives the matrix of `main_paper.tex:2115` a unit eigenvector
`z_ℓ` with eigenvalue `1`. Route: `1 - secMat` is singular, so its kernel is nonzero
(`Matrix.exists_mulVec_eq_zero_iff`); normalize any kernel vector. -/
theorem exists_unit_eigvec_secMat {M r : ℕ} {rk : Fin M → ℕ}
    {θ : (i : Fin M) → Fin (rk i) → ℝ}
    {R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ} {w : Fin M → ℝ} {γ : ℝ}
    (h : IsSecularRoot θ R w γ) :
    ∃ z : EuclideanSpace ℝ (Fin r), ‖z‖ = 1 ∧
      secMat θ R w γ *ᵥ WithLp.ofLp z = WithLp.ofLp z

/-- One summand of the limit of `main_paper.tex:2119`. -/
noncomputable def swTerm {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w c : Fin M → ℝ) (γ : ℝ)
    (z : EuclideanSpace ℝ (Fin r)) : ℝ :=
  (1 - ∑ i, c i * w i ^ 4 / (γ - w i ^ 2) ^ 2) /
    (γ * (WithLp.ofLp z ⬝ᵥ (secDerivMat θ R w γ *ᵥ WithLp.ofLp z)))

/-- The right side of `main_paper.tex:2119`, the limit of
`prop:gen_rank_stacksvd_singleweight`. -/
noncomputable def swLimit {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w c : Fin M → ℝ)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r)) : ℝ :=
  ∑ l, swTerm θ R w c (γ l) (z l)

structure EigSep {M r : ℕ} {rk : Fin M → ℕ}
    (θ : (i : Fin M) → Fin (rk i) → ℝ)
    (R : (i : Fin M) → Matrix (Fin r) (Fin (rk i)) ℝ) (w c : Fin M → ℝ)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r)) : Prop where
  /-- each `γ_ℓ` is a root of the secular equation above `sup_i w_i²` -/
  root : ∀ l, IsSecularRoot θ R w (γ l)
  /-- the roots are listed in strictly decreasing order, so no two of them are tied -/
  sorted : StrictAnti γ
  /-- `z_ℓ` is a unit eigenvector of the matrix of `main_paper.tex:2115` at eigenvalue `1` -/
  eigvec : ∀ l, ‖z l‖ = 1 ∧ secMat θ R w (γ l) *ᵥ WithLp.ofLp (z l) = WithLp.ofLp (z l)
  /-- the detectability threshold of `main_paper.tex:2106` -/
  thresh : ∀ l, ∑ i, c i * w i ^ 4 / (γ l - w i ^ 2) ^ 2 < 1
```

The performance functional (`RankR/SingleWeight/Defs.lean`):

```lean
/-- The performance of single-weight stacksvd (`main_paper.tex:2089`): `‖V̂_stacksvd(w)ᵀ V‖_F²`
in index projector form, `∑_l ∑_k overlapIdx (X_stack(w)) l v_k`. At a simple `l`-th
eigenvalue each summand is `(v̂_lᵀ v_k)²` (`perfSW_eq_inner_sq`). -/
noncomputable def perfSW (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (N : ℕ)
    (ω : Ω N) : ℝ :=
  ∑ l : Fin r, ∑ k : Fin r, overlapIdx (m.stackXW w N ω) (l : ℕ) (m.colVecG N k)

/-- `v̂_l`, the `l`-th right singular vector of the weighted stack `X_stack(w)`, with the sign
fixed by the paper's convention `⟪v̂_l, v_l⟫ ≥ 0` exactly as `vhatStackR` fixes it
(`RankR/StackGamma.lean:364`). -/
noncomputable def vhatSW (m : UnalignedModelR μ M n d r rk) (w : Fin M → ℝ) (l : Fin r)
    (N : ℕ) (ω : Ω N) : EuclideanSpace ℝ (Fin (d N)) :=
  if 0 ≤ ⟪vEig (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) (l : ℕ),
      m.colVecG N l⟫_ℝ then
    vEig (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) (l : ℕ)
  else
    -vEig (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) (l : ℕ)
```

Layer 1 (`RankR/SingleWeight/Main.lean`):

```lean
structure SingleWeightLaw (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r)) : Prop where
  /-- per pair `(l, k)`: the squared overlap of `v_k` with the `l`-th right singular subspace
  of `X_stack(w)` tends to `swTerm … (γ l) (z l) * (z l)_k²` -/
  align : ∀ l k : Fin r, TendstoInProb μ
    (fun N ω => overlapIdx (m.stackXW w N ω) (l : ℕ) (m.colVecG N k))
    (SingleWeight.swTerm (fun i => (m.tbl i).θ) m.R w c (γ l) (z l) *
      (WithLp.ofLp (z l) k) ^ 2)
  /-- the `l`-th eigenvalue of the weighted stack Gram matrix is almost surely simple at that
  index -/
  simpleIdx : ∀ (l : Fin r) (N : ℕ), ∀ᵐ ω ∂(μ N),
    SimpleIdx (m.stackGramW w N ω) (m.isHermitian_stackGramW w N ω) (l : ℕ)

theorem prop_gen_rank_stacksvd_singleweight
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r))
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z)
    (law : m.SingleWeightLaw w c γ z) :
    TendstoInProb μ (fun N ω => m.perfSW w N ω)
      (SingleWeight.swLimit (fun i => (m.tbl i).θ) m.R w c γ z)

theorem prop_gen_rank_stacksvd_singleweight_inner
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r))
    (law : m.SingleWeightLaw w c γ z) (l k : Fin r) :
    TendstoInProb μ (fun N ω => ⟪m.vhatSW w l N ω, m.colVecG N k⟫_ℝ ^ 2)
      (SingleWeight.swTerm (fun i => (m.tbl i).θ) m.R w c (γ l) (z l) *
        (WithLp.ofLp (z l) k) ^ 2)
```

The Gaussian discharge and the two facades (`RankR/SingleWeight/Het/Sup.lean`):

```lean
theorem singleWeightLaw_of_gaussian [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (hrn : ∀ N, r ≤ ∑ i with w i ≠ 0, n i N)
    {γ : Fin r → ℝ} {z : Fin r → EuclideanSpace ℝ (Fin r)}
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z) :
    m.SingleWeightLaw w c γ z

theorem singleWeightLaw_of_gaussian_shift [NeZero M] [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    {γ : Fin r → ℝ} {z : Fin r → EuclideanSpace ℝ (Fin r)}
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z) :
    ∃ k, (m.shift k).SingleWeightLaw w c γ z

theorem prop_gen_rank_stacksvd_singleweight_gaussian [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r))
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z) :
    TendstoInProb μ (fun N ω => m.perfSW w N ω)
      (SingleWeight.swLimit (fun i => (m.tbl i).θ) m.R w c γ z)

theorem prop_gen_rank_stacksvd_singleweight_inner_gaussian [NeZero M]
    [∀ N, IsProbabilityMeasure (μ N)]
    (m : UnalignedModelR μ M n d r rk) (w c : Fin M → ℝ) (hc : ∀ i, 0 < c i)
    (hreg : ∀ i, (m.tbl i).Regime (c i)) (hG : m.JointGaussianNoise)
    (γ : Fin r → ℝ) (z : Fin r → EuclideanSpace ℝ (Fin r))
    (hsep : SingleWeight.EigSep (fun i => (m.tbl i).θ) m.R w c γ z) (l k : Fin r) :
    TendstoInProb μ (fun N ω => ⟪m.vhatSW w l N ω, m.colVecG N k⟫_ℝ ^ 2)
      (SingleWeight.swTerm (fun i => (m.tbl i).θ) m.R w c (γ l) (z l) *
        (WithLp.ofLp (z l) k) ^ 2)
```

The rank-one reduction that checks the new layer against the proved tree
(`RankR/SingleWeight/ScalarsOne.lean`, not quoted here): at `r = 1`, `r_i = 1`, `R_i = 1`,
`isSecularRoot_one_iff` shows `IsSecularRoot` is `Scalars.IsGammaTop` of `lem:secular_equation`
(section 4.3), and `swLimit_one` shows `swLimit` is `Scalars.Lw θ c w` of `thm:stacksvd_weighted`
(section 4.1), given `Scalars.Assumption4` (the shape difference is modeling choice 4 below).
Numeric check: the two sides agree to 3.6e-15 over the scan of
`notes/archive/singleweight_plan.md` section 2.1.

**Hypotheses in words.**

- `sigMat`, `secMat`, `secDerivMat` carry `θ` and `R` as explicit functions, the paper's
  `Θ_i` and `R_i`, at general `rk : Fin M → ℕ` (the paper's `r_i`).
- `IsSecularRoot θ R w γ`: `γ` is above `max_i w_i²` and solves the determinant equation of
  `main_paper.tex:2115`.
- `exists_unit_eigvec_secMat`: given a root `h : IsSecularRoot θ R w γ`, a unit `z` with
  `secMat θ R w γ *ᵥ z = z` exists. This is the first sentence of the proposition, and it
  carries no probability and no `SingleWeightLaw`.
- `swTerm`, `swLimit`: the summand and the sum of `main_paper.tex:2119`, read at any `γ` and
  `z`, not only at a root; `EigSep` and `SingleWeightLaw` are what make them the right
  numbers.
- `EigSep θ R w c γ z`: `root` says each `γ l` solves the secular equation; `sorted` says the
  `r` values are strictly decreasing, so index `l` is unambiguous; `eigvec` says `z l` is the
  unit eigenvector `exists_unit_eigvec_secMat` gives at `γ l`; `thresh` is the paper's
  detectability sum of `assum:gen_rank_stacksvd_eig_sep` (`main_paper.tex:2106`).
- `m : UnalignedModelR μ M n d r rk`: the model of `assum:unaligned`, general `R_i`, without
  the paper's `Rank(∑ R_i R_iᵀ) = r` (section 1).
- `w c : Fin M → ℝ`: the one weight per table and the aspect ratios.
- `perfSW m w N ω`: `‖V̂_stacksvd(w)ᵀ V‖_F²` in projector form; `vhatSW m w l N ω` is `v̂_l`,
  signed by `⟪v̂_l, v_l⟫ ≥ 0`.
- `SingleWeightLaw m w c γ z`: `align` is the per-pair random matrix theory limit of
  `main_paper.tex:2153`; `simpleIdx` gives simplicity at every sorted index `l < r`, which
  turns the projector sum into the paper's `(v̂_lᵀ v_k)²`
  (`prop_gen_rank_stacksvd_singleweight_inner`).
- `hsep : SingleWeight.EigSep …`: `assum:gen_rank_stacksvd_eig_sep`, on the theorem's own `γ`
  and `z`.
- `law : m.SingleWeightLaw w c γ z`: the random matrix theory input, a hypothesis at Layer 1
  as everywhere in this development (`CLAUDE.md` hard rule 2), and discharged for Gaussian
  noise by `singleWeightLaw_of_gaussian`.
- `hc : ∀ i, 0 < c i` and `hreg`: `eq:RMT_limit`, the regime `n_i/d → c_i`; `hG` is
  `assum:general_noise` in the Gaussian case (both as in `thm_rank_r_stacksvd_gaussian`,
  section 6.7).
- `hrn : ∀ N, r ≤ ∑ i with w i ≠ 0, n i N` (on `singleWeightLaw_of_gaussian` only): at
  every `N` the tables of nonzero weight have at least `r` rows in total, so that the
  `d N × d N` Gram matrix of the weighted stack can have `r` simple nonzero eigenvalues;
  without it the `simpleIdx` field is false at an `N` with fewer rows. The two facades take
  no such hypothesis: they are limits, proved on the shifted model `m.shift k` (where the
  regime gives the bound at every `N`) and carried back by
  `SpikedModel.tendstoInProb_of_shift` (`notes/FLAGGED.md` D38).
- The weights carry no hypothesis, as in the paper; see below for how the Gaussian
  discharge handles a zero weight.

**Hypotheses beyond the paper.** `EigSep.sorted : StrictAnti γ` strengthens the paper's "the
top `r` eigenvalues of `G` … are distinct" (`assum:gen_rank_stacksvd_eig_sep`) to a strict
order on the sequence `γ`, not only pairwise inequality; this is finding E11 (section 7.2):
without it the limit is not determined by `(γ, z)` alone, because at a tied root the
performance (a sum of reciprocals) changes with which basis of the tied eigenspace is used
for `z`, unless that basis is chosen orthogonal for the quadratic form `secDerivMat`
(`notes/archive/singleweight_plan.md` scan 3.1). `sorted` is restrictive in that narrow sense (it
excludes the tie-tolerant basis condition the plan calls `korth`, not stated here) and
derivable from the paper's own clause otherwise. `EigSep.eigvec` is paper-implicit: the
paper's `z_ℓ` is "a unit eigenvector"; making it an explicit field of the hypothesis, rather
than reconstructing it inside the theorem, is a naming choice, not a strengthening, since
`exists_unit_eigvec_secMat` proves one exists at every root. `SingleWeightLaw` itself is the
Layer 1 black box, paper-implicit from `eq:RMT_limit` and the appendix proof at
`main_paper.tex:2123` to `:2167` (Theorem 1 of `liu2023asymptotic`, applied as in
`thm:stacksvd_weighted`); no field of it adds a condition beyond that citation, and the
Gaussian discharge proves it from the model. The weights carry no condition, as in the
paper. The almost sure simplicity of the weighted Gram matrix
(`simpleSpec_ae_stackGramW`, `RankR/Het/Simplicity.lean`) needs every weight nonzero, since a
zero weight deletes a row block of the noise; the drop step of
`RankR/SingleWeight/Het/Sub.lean` handles this: the tables of weight zero contribute nothing
to the Gram matrix, so it equals the Gram matrix of the sub-model on the tables of nonzero
weight (`stackGramW_subRk`), where every weight is nonzero
(`simpleSpec_ae_stackGramW_of_exists`). The one table of nonzero weight that the alignment
limit (`align_sw_of_gaussian`) and the shift need comes from `EigSep.root` itself: at the
zero weight the secular matrix is `0` and `det (1 - 0) = 1 ≠ 0`
(`SingleWeight.EigSep.exists_ne_zero`); at `r = 0` every field of the law is vacuous.

**Modeling choices** (`notes/archive/singleweight_plan.md` sections 4.1 and 4.2).

1. `perfSW` is the index-projector sum, not the paper's `r`-frame quantity. A top-`r` frame of
   the weighted stack's eigenspace does not exist at a tie, so the projector sum is the total
   object; it equals the paper's `‖V̂_stacksvd(w)ᵀ V‖_F²` once each of the `r` eigenvalues is
   simple at its own sorted index (`perfSW_eq_inner_sq`), which `SingleWeightLaw.simpleIdx`
   supplies almost surely. Same choice as `frobSqStackR` (`RankR/StackGamma.lean:422`).
2. `γ` and `z` are arguments of the theorem, not definitions computed from the model. Defining
   them would need a choice for `z` at a tie and a root-counting argument for `γ`; passing
   them through `EigSep` keeps the statement total and choice-free, and it is what makes the
   distinctness clause of scan 3.1 unnecessary as a separate hypothesis. `EigSep.sorted` fixes
   the correspondence between index `l` and the sorted position of the outlier.
3. The eigenvector half of the proposition is the separate deterministic theorem
   `exists_unit_eigvec_secMat`. It carries no probability, so it is proved and used on its
   own rather than folded into the limit theorem.
4. The threshold clause `thresh` lives in `EigSep`, not inside `swLimit`. The alternative
   (matching `Scalars.Lw`, which is `0` below `eq:assumption4`) is closer to the truth away
   from the hypothesis (plan scan 3.2) but further from the paper's own statement; this is
   why `swLimit_one` needs the side condition `Scalars.Assumption4` to match `Lw`.
5. `SingleWeightLaw.align` is stated in **projector** form (`overlapIdx`), which needs no
   simplicity and is stronger than the paper's inner-product form by Bessel
   (`overlapIdx_ge_inner_sq`); `simpleIdx` is what recovers the inner form
   (`prop_gen_rank_stacksvd_singleweight_inner`). Same shape as `HeteroLawR.align`
   (`RankR/StackGamma.lean:473`, section 2.5).

**Proof sketch.** `exists_unit_eigvec_secMat` (`Scalars.lean:94`) is the paper's own route:
`1 - secMat θ R w γ` is singular by `h.2`, so `Matrix.exists_mulVec_eq_zero_iff` gives a
nonzero kernel vector, and normalizing it gives the unit eigenvector at eigenvalue `1`.
`prop_gen_rank_stacksvd_singleweight` (`Main.lean:89`) is three lines: sum `law.align` over
the `r²` pairs `(l, k)` with `TendstoInProb.finsum` (`Prob/TendstoInProb.lean:396`); for each
`l` the inner sum over `k` of `(z l)_k²` is `‖z l‖² = 1` by `hsep.eigvec`, so the double sum
collapses to `∑_l swTerm … (γ l) (z l) = swLimit …`. Same shape as
`thm_rank_r_stacksvd_frobenius` (`RankR/StackMain.lean:165`). The inner form
(`Main.lean:121`) reads `law.align l k` on the almost sure event `law.simpleIdx l N`, through
`normSq_specProjIdx_eq_inner_sq`; the necessity scan found it needs no `EigSep`
binder, since the pairwise limit already holds under the law alone (its mirror
`thm_rank_r_stacksvd_inner`, `RankR/StackMain.lean:150`, carries none either).

The Gaussian discharge (`RankR/SingleWeight/Het/`, 8 files, 4387 lines, Track G) mirrors
the chain of `RankR/Het/` (section 6.7) at general `R_i` and replaces the per-column
secular scalars by the matrix ones: (1) `Het/Scalars.lean` defines the matrix `F(x)` of
resolvent limits (`swFmat`), the outlier `ρ_l = swRho c w (γ l)` of a root `γ_l` (the
`MPhet.zfun` image of `-1/γ_l`), and `ν_l = swNu …` (the reciprocal overlap), and proves
`1/(ρ_l ν_l) = swTerm … (γ_l) (z_l)`, the paper's summand; `bHet_lt_swRho` puts every root's
outlier above the bulk edge from the paper's threshold clause. (2) `Het/Forms.lean` proves
the resolvent limits of the columns of the stacked signal matrix (`ResolventLimitsSW`,
`resolventLimitsSW_of_gaussian`) and their Gram limit at general `R_i`, from the
heteroscedastic resolvent theory of `RankR/Het/Forms.lean`. (3) `Het/Frame.lean` and
`Het/Count.lean` are deterministic: a frame `Z` of approximate eigenvectors with a small
residual and a Gram matrix near the identity locates one eigenvalue near each outlier
(`residual_bound_Z`, `gram_bound_Z`), and the count of eigenvalues above a gap point equals
the number of outliers above it (`count_of_split_of_frame_gap`), which also reads the
overlap of a column off the frame (`align_detZ`). (4) `Het/Outliers.lean` applies the frame
`Z = zMat z` of the eigenvectors of `EigSep` to the sample matrix: the count above `τ`
(`tendsto_measure_count_Ioi_tau_sw`) and the projected norm of column `k`
(`tendstoInProb_normSq_specProj_Ioi_tau_sw`). (5) `Het/Align.lean` reads the sorted
eigenvalue `l` (`tendstoInProb_eigVal_sw`) and the projector overlap at index `l`
(`align_sw_of_gaussian`) from two counts at `ρ_l ∓ δ`, as `RankR/Het/Align.lean` does.
(6) `Het/Sub.lean` (F18c) drops the tables of weight zero: the sub-model `subRk` on the
support of `w` has the same weighted Gram matrix (`stackGramW_subRk`), so the almost sure
simplicity of `RankR/Het/Simplicity.lean`, which needs every weight nonzero, transfers
(`simpleSpec_ae_stackGramW_of_exists`); and `EigSep.exists_ne_zero` supplies the one table
of nonzero weight from the secular equation. (7) `Het/Sup.lean` assembles the structure:
`align` from (5), `simpleIdx` from (6), the side condition `d N = p N + r` from the tail
shift `exists_shift_sw`, and the facades from Layer 1 on the shifted model.

**Scope.** Layer 1 and the Gaussian discharge. The eigenvector half
(`exists_unit_eigvec_secMat`) needs no discharge: it is deterministic. `SingleWeightLaw`
has an explicit witness through the discharge for any Gaussian model that meets `EigSep`;
no separate rank-one instance is written (section 8.1): the natural one needs a bridge from
`MultiTableModel` to `UnalignedModelR` that does not exist yet (F18a's bridge item). The
discharge alone does not close `prop:singleweight_suboptimality` (section 6.9): its `hlaw`
asks for the limit at every positive weight pair of the witness, where `EigSep` fails at
extreme weight ratios and the paper's `swLimitEx` has undetectable branches (risk 4 of
`notes/archive/trackG_plan.md`). Item F18b (2026-09-07) proves those branches separately
(`RankR/SingleWeight/Tie.lean`, `Regimes.lean`, `Het/OneOutCount.lean`, `Het/OneOutDet.lean`,
`Het/OneOutAlign.lean`, `Het/OneOutBulk.lean`, `Suboptimality.lean`), so `hlaw` is a theorem
(`hlaw_witness`) and section 6.9 is unconditional.

**Numeric check.** `scripts/numeric/check_singleweight.py` (seed 20260905, numpy only): the finite-`d`
deterministic identity `det(G - γI) = det(Σ - γI) det(I_r - secMat θ R w γ)` holds to between
8.9e-16 and 7.1e-15 on three instances and to 3.9e-15 worst case on 300 random instances; the
performance formula matches Monte Carlo simulation within about 1 standard error at `d = 1600`
on two general-`R_i` instances and on the paper's own suboptimality instance
(`notes/archive/singleweight_plan.md` section 2).

**Audit trail.** `notes/archive/singleweight_plan.md` (status user OK 2026-09-05, statements as
compiled), hypothesis necessity scan of the same note section 3 (all four clauses of
`assum:gen_rank_stacksvd_eig_sep` needed for the value, not only the proof route; the
distinctness finding is E11), adversarial risk review section 7 (four risks, none overturning
the statement). Proved 2026-09-05; delegation record in the plan note section 6. Gaussian
discharge: `notes/archive/trackG_plan.md` (the route, the unit table, the risks, and section 8, the
execution record with the model and the coordinator's check per unit), the numeric checks
of `notes/archive/trackG_specs/` (`check_trackG.py`, `check_trackG_mc.py`, `check_frame.py`,
`audit_G0.py`; seed 20260905), decisions D37 (`hw`, resolved by F18c the same day) and D38
(`hrn`, the shift) of `notes/FLAGGED.md`. The user OK'd the Layer 1 statements; the
Gaussian signatures were fixed by Claude (D37, D38) and wait for the user's review at the
handoff.

### 6.9 `prop:singleweight_suboptimality`

**Paper** (`main_paper.tex:912` to `:917`, the statement; `:2171` to `:2194`, the instance of
the proof).

```latex
\paragraph*{General (single) weighting} Alternatively, we can consider using a single weight per matrix. However, while we can compute the asymptotic performance as a function of the weighting vector and the rotation matrices, it does not yield a closed form optimization and requires knowledge of the $R_i$ (\Cref{prop:gen_rank_stacksvd_singleweight} in \Cref{sec:insufficiency_of_global_weights_stacksvd}).
Nevertheless, through a more careful analysis of the example in \Cref{fig:unaligned_sim} (\Cref{eq:psi_equation} with $\sin(\psi)=0$), we are able to prove the following Proposition (details in \Cref{sec:appendix_insufficiency_single_weights_stacksvd}):

\begin{prop} \label{prop:singleweight_suboptimality}
    In the general rank-$r$ setting, there exists a problem instance such that unweighted \svdstack outperforms optimally weighted \stacksvd with a single weight per table.
\end{prop}
```

The instance, from the proof (`main_paper.tex:2171` to `:2194`):

```latex
We now apply the result to demonstrate that unweighted \svdstack can outperform optimally weighted \stacksvd with a single weight per table in the general rank setting. Take 
\begin{equation}
    R_1 = \begin{bmatrix}
        1 \\ 0
    \end{bmatrix}, \quad R_2 = \begin{bmatrix}
        0 \\ 1
    \end{bmatrix}, \Theta_1 = \Theta_2 = \diag(\theta_0), \quad \theta_0^4 > c_0, \quad c_1 = c_2 = c_0.
\end{equation}
To derive the performance of optimally weighted \svdstack (which coincides with unweighted \svdstack in this setting), we may apply Theorem \ref{thm:gen_rank_weight_svdstak} with
\begin{equation*}
    A_{\beta, R} = \begin{bmatrix}
        1 & 0 \\ 0 & 1
    \end{bmatrix}, \quad (W\opt)^{-2} = \begin{bmatrix}
        1 - \beta_0^2 & 0 \\ 0 & 1 - \beta_0^2
    \end{bmatrix}
\end{equation*}
and therefore $\|\Vsvdstack(W\opt)^\top V\|_F^2 \ip 2 - (1 - \beta_0^2) - (1 - \beta_0^2) = 2 \beta_0^2$. For \stacksvd, we first assume without loss of generality that $w_1 > w_2$ and that assumption \ref{assum:gen_rank_stacksvd_eig_sep} holds. Otherwise, at least one component is undetectable and the performance is bounded above by $\beta_0^2$. In this case, a simple application of \Cref{prop:gen_rank_stacksvd_singleweight} yields
\begin{equation*}
\begin{split}
        \gamma_{\ell} &= w_{\ell}^2(1+ \theta_0^2), \quad w_2^2(1+\theta_0^2) > w_1^2, \quad z_{\ell} = e_{\ell} \in \R^2 \\
        \| \Vstacksvd(w_1, w_2)^{\top} V\|_F^2 &\ip \frac{\theta_0^4 - c_0 - \frac{c_0 w_2^4 \theta_0^4}{(w_1^2(1 + \theta_0^2) - w_2^2)^2}}{\theta_0^2(1 + \theta_0^2)} + \frac{\theta_0^4 - c_0 - \frac{c_0 w_1^4 \theta_0^4}{(w_2^2(1 + \theta_0^2) - w_1^2)^2}}{\theta_0^2(1 + \theta_0^2)} < 2 \beta_0^2.
\end{split}
\end{equation*}
When $w_1 = w_2$, the performance is obtained by a result for unweighted \stacksvd, and is also strictly less than $2 \beta_0^2$. This proves \Cref{prop:singleweight_suboptimality}, showing that in the general rank-$r$ setting, there exists a problem instance such that unweighted \svdstack outperforms optimally weighted \stacksvd with a single weight per table.
```

**Lean.** Two forms: `prop_singleweight_suboptimality_of_law`, with the hypothesis `hlaw`
(the single-weight limit on the witness at every positive weight pair) in place of its
proof, and `prop_singleweight_suboptimality_gaussian`, unconditional, which discharges
`hlaw` by `hlaw_witness` (`RankR/SingleWeight/Suboptimality.lean`, F18b, 2026-09-07). The witness
model is a Gaussian model (`RankR/SingleWeight/Example.lean`,
`RankR/SingleWeight/Existence.lean`). The instance's scalars:

```lean
/-- `γ_ℓ = w_ℓ²(1 + θ_0²)`, the `ℓ`-th secular root of the instance
(`main_paper.tex:2190`). -/
noncomputable def gammaEx (θ₀ : ℝ) (w : Fin 2 → ℝ) (l : Fin 2) : ℝ := w l ^ 2 * (1 + θ₀ ^ 2)

noncomputable def swTermEx (θ₀ c₀ : ℝ) (w : Fin 2 → ℝ) (l : Fin 2) : ℝ :=
  (θ₀ ^ 4 - c₀ - c₀ * w (other l) ^ 4 * θ₀ ^ 4 /
      (gammaEx θ₀ w l - w (other l) ^ 2) ^ 2) / (θ₀ ^ 2 * (1 + θ₀ ^ 2))

/-- Root `ℓ` of the instance is detectable: `γ_ℓ = w_ℓ²(1 + θ_0²)` lies above every `w_i²`,
and the threshold sum of `assum:gen_rank_stacksvd_eig_sep` (`main_paper.tex:2106`) is below
`1`. The two clauses are `IsSecularRoot`'s first clause and `EigSep.thresh` at `c = c_0`,
written on the instance. -/
def DetectableEx (θ₀ c₀ : ℝ) (w : Fin 2 → ℝ) (l : Fin 2) : Prop :=
  Scalars.wSqMax w < gammaEx θ₀ w l ∧
    ∑ i, c₀ * w i ^ 4 / (gammaEx θ₀ w l - w i ^ 2) ^ 2 < 1

open Classical in
/-- The closed form of the instance, total in `w` (`main_paper.tex:2190` to `:2191`). Three
branches, in this order.

1. Both roots detectable: the sum of the two terms of the paper's display.
2. Exactly one root detectable: that one term. The other component sits in the bulk and
   contributes nothing (`main_paper.tex:2187`).
3. Neither root detectable: `0`.

The `if` conditions are `DetectableEx θ₀ c₀ w 0` and `DetectableEx θ₀ c₀ w 1`, in that
order. -/
noncomputable def swLimitEx (θ₀ c₀ : ℝ) (w : Fin 2 → ℝ) : ℝ :=
  if DetectableEx θ₀ c₀ w 0 then
    if DetectableEx θ₀ c₀ w 1 then
      swTermEx θ₀ c₀ w 0 + swTermEx θ₀ c₀ w 1
    else
      swTermEx θ₀ c₀ w 0
  else
    if DetectableEx θ₀ c₀ w 1 then
      swTermEx θ₀ c₀ w 1
    else
      0

theorem swLimitEx_lt : ∀ w : Fin 2 → ℝ, 0 < w 0 → 0 < w 1 →
    swLimitEx (8 / 5) 1 w < 2 * betaSq (8 / 5) 1
```

`other : Fin 2 → Fin 2` is the other index (`other 0 = 1`, `other 1 = 0`), and
`swLimitEx_eq_swLimit` (`Example.lean:200`, not quoted here) proves that on the two-root
branch this closed form is the general `swLimit` of section 6.8 read at `γ_ℓ`, `z_ℓ = e_ℓ`,
`R = Rone (RankR.Example.Rex 0)`: this is the step that lets a Gaussian discharge of
`SingleWeightLaw` on the instance recover the paper's own display.

The witness model (`RankR/SingleWeight/Existence.lean`): `M = 2` Gaussian tables, `r = 2`,
one spike each, `θ_0 = 8/5`, `R = Rone (RankR.Example.Rex 0)` (that is `R_1 = (1,0)ᵀ`,
`R_2 = (0,1)ᵀ`), `c_1 = c_2 = 1` (row counts `n_i N = d N = N + 3`). Its svdstack half is
already unconditional:

```lean
theorem witness_perfRG_tendsto :
    TendstoInProb mu (fun N ω => mdl.perfRG N ω) (2 * betaSq (8 / 5) 1)
```

and the existence proposition itself:

```lean
theorem prop_singleweight_suboptimality_of_law
    (hlaw : ∀ w : Fin 2 → ℝ, 0 < w 0 → 0 < w 1 →
      TendstoInProb mu (fun N ω => mdl.perfSW w N ω) (swLimitEx (8 / 5) 1 w)) :
    ∃ (Ω : ℕ → Type) (_ : ∀ N, MeasurableSpace (Ω N)) (μ : ∀ N, Measure (Ω N))
      (_ : ∀ N, IsProbabilityMeasure (μ N)) (n : Fin 2 → ℕ → ℕ) (d : ℕ → ℕ)
      (m : UnalignedModelR μ 2 n d 2 (fun _ => 1)),
      (∀ i j, (m.tbl i).θ j = 8 / 5) ∧ (∀ i, (m.tbl i).Regime 1) ∧
      m.JointGaussianNoise ∧ m.R = Rone (RankR.Example.Rex 0) ∧
      TendstoInProb μ (fun N ω => m.perfRG N ω) (2 * betaSq (8 / 5) 1) ∧
      (∀ w : Fin 2 → ℝ, 0 < w 0 → 0 < w 1 →
        TendstoInProb μ (fun N ω => m.perfSW w N ω) (swLimitEx (8 / 5) 1 w) ∧
        swLimitEx (8 / 5) 1 w < 2 * betaSq (8 / 5) 1) ∧
      2 * betaSq (8 / 5) 2 < 2 * betaSq (8 / 5) 1
```

and the unconditional form (`RankR/SingleWeight/Suboptimality.lean`):

```lean
theorem prop_singleweight_suboptimality_gaussian :
    ∃ (Ω : ℕ → Type) (_ : ∀ N, MeasurableSpace (Ω N)) (μ : ∀ N, Measure (Ω N))
      (_ : ∀ N, IsProbabilityMeasure (μ N)) (n : Fin 2 → ℕ → ℕ) (d : ℕ → ℕ)
      (m : UnalignedModelR μ 2 n d 2 (fun _ => 1)),
      (∀ i j, (m.tbl i).θ j = 8 / 5) ∧ (∀ i, (m.tbl i).Regime 1) ∧
      m.JointGaussianNoise ∧ m.R = Rone (RankR.Example.Rex 0) ∧
      TendstoInProb μ (fun N ω => m.perfRG N ω) (2 * betaSq (8 / 5) 1) ∧
      (∀ w : Fin 2 → ℝ, 0 < w 0 → 0 < w 1 →
        TendstoInProb μ (fun N ω => m.perfSW w N ω) (swLimitEx (8 / 5) 1 w) ∧
        swLimitEx (8 / 5) 1 w < 2 * betaSq (8 / 5) 1) ∧
      2 * betaSq (8 / 5) 2 < 2 * betaSq (8 / 5) 1
```

**Hypotheses in words.**

- `θ₀ c₀ : ℝ`, `w : Fin 2 → ℝ`: the instance's shared spike strength, aspect ratio, and the
  two table weights; `8/5` and `1` below are the paper's `θ_0` and `c_0`, chosen so every
  constant is rational (`θ_0² = 64/25`, `θ_0⁴ = 4096/625 > 1`, `β_0² = 39/64`,
  `2β_0² = 39/32`).
- `gammaEx`, `swTermEx`, `DetectableEx`, `swLimitEx`: the closed form of
  `main_paper.tex:2190` to `:2191`, made total (modeling choice 1 below).
- `swLimitEx_lt`'s `0 < w 0` and `0 < w 1`: both table weights are positive, the paper's
  implicit reading of "a single weight per table" (modeling choice 2 below, and see
  "Hypotheses beyond the paper").
- `mdl : UnalignedModelR mu 2 nn dd 2 (fun _ => 1)`: the witness model, a real Gaussian model,
  not a hypothesis.
- `witness_perfRG_tendsto`: no hypothesis beyond the model. The svdstack half of the claim is
  unconditional.
- `hlaw`: the stacksvd limit on the instance, for every positive weighting, the one
  hypothesis of the `_of_law` form. It is what the paper's appendix proof reads off
  `prop_gen_rank_stacksvd_singleweight` on `mdl`. `hlaw_witness` proves it in three
  regimes. Both secular roots detectable and `w 0 ≠ w 1`: the Gaussian discharge of
  `SingleWeightLaw` (section 6.8, Track G) applies, since `Regimes.lean` shows `EigSep` on
  the witness there (`eigSep_of_detectable`) and `swLimit = swLimitEx` (`hlaw_both`).
  Exactly one root detectable (an extreme weight ratio, where `EigSep` fails and
  `swLimitEx` has one term): the count of eigenvalues above the edge is at most one
  (`Het/OneOutCount.lean`), the one outlier converges and its eigenvector aligns with the
  detectable direction (`Het/OneOutAlign.lean`), and the top-2 window at the edge carries no
  mass of the spike directions (`Het/OneOutBulk.lean`), so the second index contributes `0`
  (`align_bulk_sw_one`) and the sum is `swTermEx l₀` (`hlaw_one`). The tie `w 0 = w 1`: the
  weighted stack is a multiple of the unweighted stack, and
  `prop_stacksvd_subspace_general_gaussian` gives `swLimitEx`, as the paper's last sentence
  says (`hlaw_tie`, `Tie.lean`, `swLimitEx_tie`).
- The conclusion existentially quantifies the whole model, in the shape of
  `remark_svd_outperform_stack_exists` (section 5.6): a reader gets a concrete `Ω`, `μ`, `n`,
  `d`, `m` with no hypothesis left over (the `_gaussian` form; `hlaw` alone in the `_of_law`
  form).

**Hypotheses beyond the paper.** `0 < w 0` and `0 < w 1`, in `swLimitEx_lt` and in the
existence conclusion, are restrictive against the paper's own phrase "a single weight per
table", which places no positivity condition on `w`. `w_i = 0` is fine for the closed-form
value `swLimitEx` on its own (`notes/archive/singleweight_plan.md` scan 3.6), but on this instance it
drops table `i` entirely from the weighted stack, and with `R_1 = (1,0)ᵀ`, `R_2 = (0,1)ᵀ` the
surviving table spans only one direction, so the secular equation has only 1 root where `r =
2` are needed for the sum `swLimit` to be defined at all (plan scan 3.6, the same failure the
rank condition of `assum:unaligned` triggers at scan 3.4). That is why the existence statement
requires `0 < w i` (plan section 4.3, modeling choice 5). `hlaw` of the `_of_law` form is
paper-implicit: it is exactly what the appendix proof at `main_paper.tex:2171` to `:2194`
establishes, once `prop:gen_rank_stacksvd_singleweight` is available, and `hlaw_witness`
proves it (F18b). One difference from the paper's text: where exactly one root is
detectable the paper states only the bound `β_0²` ("bounded above by", `:2187`), and the
Lean proves the exact limit `swTermEx l₀ < β_0²` (`notes/paper_edits.md` E12).

**Modeling choices** (`notes/archive/singleweight_plan.md` sections 4.3 and 7).

1. `swLimitEx` is total in `w`, with three branches (both roots detectable, one, neither),
   rather than only the paper's displayed case `w_1 > w_2` with both detectable. The plan's
   numeric check 2.4 finds the supremum of the paper's own display sits at `w_1 = w_2`, the
   point the display's "without loss of generality" excludes; the paper's last sentence
   covers that point separately, by the unweighted rank-`r` law. A statement that only
   formalized the displayed branch would miss the supremum, and the existence claim would not
   follow, so `swLimitEx` covers all three branches and `swLimitEx_lt` is proved on each
   (adversarial risk 4 of the plan).
2. `0 < w i` in the existence statement, not `w ≠ 0` or unconstrained: see "Hypotheses beyond
   the paper" above.
3. `8/5` in place of a bound variable `θ_0`: `Existence.lean` and `Example.lean` use explicit
   rationals for the same reason elsewhere in the tree (for example `Sat.lean`'s
   `![Real.sqrt 5, 4]`); `8/5` keeps every derived constant of this instance rational.
4. The witness is an `UnalignedModelR`, the general-`R_i` family, not the fully aligned
   `UnalignedModel`: `perfSW` and `EigSep` live on `UnalignedModelR`, and no map from
   `MultiTableModel` supplies one yet (`SpikedModelR.hθnn` allows `θ = 0`, and
   `SpikedModel.toRankR` maps one table; the multi-table map is still to write).
5. `hlaw` is a hypothesis on the model `mdl` at every positive `w`, not the general
   `SingleWeightLaw` structure of section 6.8. It is the minimal input the argument needs:
   the proof of `swLimitEx_lt` never reads `EigSep` or `SingleWeightLaw`, only the closed-form
   value on the instance. `hlaw_witness` discharges it on `mdl`. Its one-root helper
   `align_bulk_sw_one` takes a detectability hypothesis `hnd : ¬ DetectableEx (8/5) 1 w l₁`
   that the statement OK'd on 2026-09-06 lacked; without it the claim is false where both
   roots are detectable (decision D39 of `notes/FLAGGED.md`; the endpoint is unchanged).

**Proof sketch.** `swTermEx_lt` (`Example.lean:119`, private) is the term-by-term bound: each
summand of `swLimitEx` equals `β_0²` minus a strictly positive penalty once the other table's
weight is positive, because `β_0² = (θ_0⁴ − c_0)/(θ_0²(1+θ_0²))` shares the denominator with
`swTermEx` (`notes/archive/singleweight_plan.md` section 2.4, observation 1). `swLimitEx_lt`
(`Example.lean:138`) case-splits on the two `DetectableEx` conditions: on the two-root branch
it sums two applications of `swTermEx_lt`; on the one-root branch it applies `swTermEx_lt`
once; on the no-root branch the value is `0 < 39/32` by `norm_num`. `witness_perfRG_tendsto`
(`Existence.lean:182`) routes through the `r_i = 1` facade
`prop_general_rank_unweighted_svdstack_general_gaussian_one` of the theorem quoted in section
6.2, then the scalar chain `limitRG_one_eq_limitR`, `limitR_example_beta`, `svdstackEx_zero_s`
of `RankR/Example.lean` (section 6.5), since `A_{β,R} = I_2` on this instance, the paper's own
remark. `prop_singleweight_suboptimality_of_law` (`Existence.lean:213`) assembles the model
data (`mdl_theta`, `mdl_regime`, `mdl_joint`, `mdl_R`), `witness_perfRG_tendsto`, `hlaw` and
`swLimitEx_lt` into the existential witness by `refine`, then closes the last conjunct
`2 betaSq (8/5) 2 < 2 betaSq (8/5) 1` by `norm_num [betaSq]` directly on the concrete
rationals. `hlaw_witness` (`Suboptimality.lean`) splits on `lt_trichotomy (w 0) (w 1)`: the
tie goes to `hlaw_tie`; otherwise the smaller weight's root `l₁` is detectable or not
(`DetectableEx` is decidable classically), which sends the pair to `hlaw_both` or `hlaw_one`
(the larger weight's root is detectable whenever the smaller one is, `detectableEx_of_le`,
and always is on this instance, `Regimes.lean`). `hlaw_one` writes `perfSW` as the sum over
the two sorted indices (`Fin.sum_univ_two`); index `0` follows Track G's sandwich of the
outlier between `swRho ± δ` with the count of `Het/OneOutCount.lean` in place of `EigSep`,
and index `1` is `align_bulk_sw_one`: the second eigenvalue sits in the edge window
(`Frame.eigenvalues₀_le_of_two_cols`, `R6het.lamMax_add_vecMulVec_le_of_secular_pos`, since
the undetectable root gives a positive secular value at the edge), and
`tendsto_measure_normSq_specProj_edge_gt_sw` bounds the overlap through
`overlapIdx_le_of_eigVal_ge`.

**Scope.** Both halves of the claim are theorems: the svdstack half `witness_perfRG_tendsto`
and the stacksvd half `hlaw_witness` (F18b, 2026-09-07, 7 files, 2149 lines:
`Tie.lean`, `Regimes.lean`, `Het/OneOutCount.lean`, `Het/OneOutDet.lean`,
`Het/OneOutAlign.lean`, `Het/OneOutBulk.lean`, `Suboptimality.lean`; plan
`notes/archive/F18b_plan.md`). `prop_singleweight_suboptimality_gaussian` is the
unconditional statement, with the conclusion of the `_of_law` form. The removal remark of
`main_paper.tex:901`, cited in the same paragraph as this proposition's setup, is a separate
claim and stays not formalized (question Q8, section 7.3).

**Numeric check.** `notes/archive/singleweight_plan.md` section 2.4: at `θ_0 = 1.6`, `c_0 = 1` the
general formula of section 6.8 and this closed form agree to 2.2e-16 at three values of
`w_2²`; the equal-weight tie point agrees with the unweighted rank-`r` law to 1.1e-16; Monte
Carlo at `d = 1600` is within about 1 standard error of the formula (largest deviation 0.028
at 0.9 standard errors, section 2.3). At `θ_0 = 8/5` (this instance's value) `betaSq (8/5) 1 =
39/64` and `2 betaSq (8/5) 2 = 0.999297753 < 39/32`, both exact rational checks inside the
Lean proof (`norm_num`), not floating point.

**Audit trail.** `notes/archive/singleweight_plan.md` (status user OK 2026-09-05), section 2.4 (the
closed-form numeric check), section 3.6 (the `w_i = 0` necessity scan), section 7 risk 4 (the
excluded-boundary supremum, resolved by the total `swLimitEx`). Proved 2026-09-05.

## 7. Known gaps and findings

### 7.1 What is not formalized

| Paper result | Label | Why not |
|---|---|---|
| The removal remark | `main_paper.tex:901` | not formalized, and wrong in general; see question Q8 and section 6.4. |

One more paper result was formalized at Layer 1 only until 2026-09-07: the suboptimality
corollary `prop:singleweight_suboptimality` (`main_paper.tex:915`, section 6.9), proved
under the hypothesis `hlaw` (the single-weight stacksvd limit on the two-table witness at
every positive weight pair). Item F18b proves `hlaw` (`hlaw_witness`), and
`prop_singleweight_suboptimality_gaussian` is unconditional. Its input, single-weight
rank-`r` stackSVD (`prop:gen_rank_stacksvd_singleweight`, `:2112`, the limit of stacksvd
under one weight per table, section 6.8), is proved at Layer 1 and for Gaussian noise
(`prop_gen_rank_stacksvd_singleweight_gaussian`, Track G); its Gaussian theorem needs the
separation assumption `EigSep`, which fails on the witness at an extreme weight ratio, so
F18b adds a one-component outlier law and an edge window for those ratios, and the
unweighted rank-`r` law at the tie. `SingleWeightLaw` has no separate rank-one witness in
`Sat.lean` (section 8.1); the discharge itself exhibits it on any Gaussian model that meets
`EigSep`.

Two model-level restrictions, both recorded:

1. **Noise.** Every Gaussian theorem assumes Gaussian noise. The paper's
   `assum:general_noise` asks only for two moments equal to `1` and a bounded fourth moment.
   The stackSVD Layer 1 theorems (sections 3.4, 3.6, 3.7, 4.1, 5.1 to 5.3, 6.3, 6.7) and the
   rank-`r` SVDstack theorems (6.1, 6.2, 6.4, 6.6) read the noise only through their law
   structure and, at rank `r`, through `UnalignedModelR.IndepNoise` (independent tables with
   arbitrary laws); a four-moment discharge of the same structures would give those results
   with no change to any statement. The rank-one SVDstack family (`lem_delocalization` 3.2,
   `thm_svd_stack_general` 3.5, `thm_svdstack_weighted` 4.2, the uniform bound of
   `SVDStack/Rayleigh.lean`) reads the noise through `MultiTableModel.IndepNoise`
   (independent tables with arbitrary laws, `SVDStack/Gram.lean`; L4 of the
   external statement audit, `notes/archive/L4_indepnoise.md`), so the same holds there.
   `thm_theta_est` (5.4) reads the noise through `MultiTableModel.ThetaEstLaw` (section 2.6)
   since F31 (2026-09-08): item P with a random direction (`‖E_j v̂_i‖² → c_j`) and the cross
   term (`u_jᵀ E_j v̂_i → 0`), which for Gaussian noise are Chebyshev bounds with second
   moments under `gaussianMatrix` (`noise_projection_topDir_tendsto`, `cross_term_tendsto`);
   with arbitrary marginals both limits fail, so a non-Gaussian discharge needs a moment
   hypothesis, not only independence (`notes/FLAGGED.md`, D35, superseded in its Layer 1
   half by F31). Until F31 the Layer 1 form took `hG : m.JointGaussianNoise` itself, the
   one rank-one exception (user decision Q3; external statement audit finding 8). Since
   2026-09-08 the two limits are also proved for any fixed law with mean `0`, variance `1`
   and a finite fourth moment (`thetaEstLaw_of_general`, section 5.4).
2. **Ordering of spikes inside a table.** `SpikedModelR.hθanti` and `hθnn` order the
   `θ_ij` strictly and make them nonnegative (F8: a zero last spike is allowed). The paper assumes distinct entries for svdstack
   only. See sections 6.3 and 6.7.

**The Furedi-Komlos exponent.** The four-moment edge of section 3.1 rested on one
combinatorial count: a closed walk of length `2k` on the complete bipartite graph
`Fin n × Fin d`, with no entry visited exactly once and with one vertex fewer than a tree
walk, costs at most `(2k)^12 / min(n, d)` per missing vertex (the Furedi-Komlos count;
Furedi and Komlos 1981 page 237, Anderson, Guionnet and Zeitouni Lemma 2.1.23). It was
proved 2026-09-10 evening in three new modules of `RMT/General/Edge/`: `Code.lean`,
`Dyck.lean` and `CountBound.lean`. The route codes each step of a walk and bounds the bad
ones by a credit-debt induction (`Code.lean`), bounds the tree-walk cell from below by a
Narayana double rotation (`Dyck.lean`), and assembles the two into the count
(`CountBound.lean`: `cellCard_mul_le_pos`, `cellCard_mul_le`, `clsCard_mul_le`).
`docs/SORRIES.md` has no row, and everything in the Stage 3 chain is proved. The count is
checked by exact enumeration for `k` up to 7 and for `(n, d)` from `(k+1, k+1)` to
`(10^12, 10^9)`: 0 violations, and it still holds with `(2k)^2` in place of `(2k)^12`, so the
exponent has a margin of `(2k)^10`.

Every result is stated as convergence **in probability**. The paper says the same. Nothing is
stated almost surely. Two facts inside the development are almost sure and finite-`N`, and
they are strengthenings, not gaps: `SingleTableLaw.topSimple` and `TableLawR.simple`.

**Two paper results that an auditor will grep and not find.** `prop:maxrow`
(`main_paper.tex:2628`) and `cor:domination` (`:2679`) are **commented out** in the source of
the AoS resubmission, together with the `\maxrow` estimator they describe. They are not part
of the paper's claims and nothing formalizes them.

**One paper lemma with no section here.** `lem:noise_projection_concentration`
(`main_paper.tex:2496`), `‖E_2 a‖₂² → c_2` for a unit `a` independent of `E_2`, is proved, not
assumed, as `SpikedModel.noise_projection_tendsto` in `ThetaEst.lean`. It is item P of the
roadmap and it is used only inside `thm:theta_est` (section 5.4).

**Well-formedness hypotheses that remain.** These are not restrictions on the mathematics;
they make an object exist or an index be in range. An auditor should read each one and check
that it is not doing hidden work.

| Hypothesis | Where | What it does |
|---|---|---|
| `[∀ N, IsProbabilityMeasure (μ N)]` | every Gaussian theorem | each `μ N` is a probability measure. Implicit in the paper. Where `JointGaussianNoise` is present it is derivable, and several theorems (for example `prop_dominance_gaussian`) do without it. |
| `[NeZero M]`, `hM : 0 < M` | the rank-`r` results, and the 26 rank-one declarations listed in `notes/archive/F22_nezero.md` (every one writes the binder in its quoted source text; see section 0 item 3) | at least one table, so the stack has rows and `r ≤ r̃` holds. |
| `[NeZero S.card]`, `S.Nonempty` | `cor.2`, section 5.3 | the kept subset is nonempty; the paper's "assuming at least one table is above the threshold". |
| `hr : 0 < r` | the rank-`r` theorems of `RankR/` that read a top-`r` eigenframe | the shared subspace is nontrivial, so the frame exists. At `r = 0` both sides are `0`, which is why the `r_i = 1` form of section 6.4 drops it. It is also absent from eight declarations that never read it (`topGap_gram_whp`, `topGap_gramW_whp`, `topGap_gramG_whp`, `topGap_gramWG_whp` and the four `prop_*_frobenius_eig*` that forwarded it); it stays everywhere else. Note that `hr : 0 < r` in `LinAlg/Eigen.lean`, `LinAlg/SpecProjPerturb.lean`, `RMT/MP7.lean` and `RMT/T.lean` binds a different `r` (a spectral rank of a bare matrix, or a real number), not the model's shared dimension. |
| `hrr : r ≤ rtot rk`, `hrM : r ≤ M` | sections 6.2, 6.4 | the paper's `r ≤ r̃`; the top-`r` spectral objects of `A_{β,R}` need it. |
| `hc : 0 < c i` | almost everywhere | `eq:RMT_limit`'s `c_i ∈ (0, ∞)`. Section 6.3 needs only `0 < ∑ i, c i`. |
| `hβdef`, `hθ`, `hR` | many | naming hypotheses. They tie a free variable to the model and prove nothing. |

### 7.2 The paper findings E1 to E12 (`notes/paper_edits.md`)

Each item is a proposal for the author. Every one names a Lean anchor.

| # | Result | Finding |
|---|---|---|
| E1 | `thm:svdstack_weighted` (line 502) | the optimality claim holds over **all nonzero** weight vectors, with no eigengap condition. Anchor: `thm_svdstack_weighted_gaussian_opt_full` (`SVDStack/Rayleigh.lean`). |
| E2 | `thm:gen_rank_weight_svdstak` (line 893) | `β_ij > 0` and the footnote can go; the upper bound holds for every `W` and is attained at `W⋆`, with no rank condition. Anchor: `thm_gen_rank_weight_svdstak_general_r_full_gaussian`. |
| E3 | `eq:stacksvd_gammak` (line ~2331) | `γ_j` should be defined below the threshold too (value `0`), so the corollary reads as one statement. Anchor: `Scalars.gammaR`. |
| E4 | `thm:rank_r_svdstack` (line ~2400) | it is not a specialization of `thm:gen_rank_weight_svdstak`: the component clause needs its own proof. Item 2 of E4 (asking for `S_j > 0`) should itself be revised; the clause holds at `S_j = 0` with limit `0`. Anchors: `thm_rank_r_svdstack_component_gaussian`, `thm_rank_r_svdstack_aggregate_gaussian`. |
| E5 | minor items | `thm:simple_thm1` is valid for every `M ≥ 1`; the strict clauses of `thm:stacksvd_binary_optimal_svd_stack` and `prop:dominance` need no change; the MLE appendix should name the marginalization as a standard fact; the two remarks need no change. |
| E6 | `remark:svd_outperform_stack` (line 603) | the `M = 3` example sits exactly on the threshold, so binary stackSVD reaches `62/72 = 0.861` against svdstack's `6/7 = 0.857`, the reverse of the claim. `θ₃ = (c₃+1)^{1/4}` restores it. Anchor: `remark_svd_outperform_stack_three_binary`. |
| E7 | `thm:rank_r_stacksvd` (line ~2306), `alg:rank_r_stacksvd` (line 2382) | the index `ℓ_j` is wrong at unordered `θ`; an instance with closed-form outlier locations and Monte Carlo limits contradicts the first display there (numerical evidence, not a certified refutation). Inside the ordered class it is right. Anchors: `Scalars.ellR`, `Scalars.ellSup`, `Scalars.ellSup_eq_ellR`. |
| E8 | `thm:stacksvd_weighted` (line 463), `prop:dominance` (line 634) | both are false at `θ ≡ 0`: the stated optimal weights are all zero, the weighted stack is the zero matrix, and `(vᵀ v̂)²` is `1`, not the claimed `γ⋆ = 0`. Anchor: `hθ : ∃ i, θ_i ≠ 0` in `thm_stacksvd_weighted_gaussian` and `prop_dominance_gaussian`; `notes/REVISION_LIST.md` R5. |
| E9 | eight errata (lines 376, 442, 1331, 1612, 842, 909, 1345, 212) | `M ≥ 2` missing; a maximum over a set that contains `∅`; a weight condition that describes an empty set; a reversed inequality; `V̂_stacksvd` not unique at a tie; two counts `2^M` that include the empty weighting; a caption that cites the wrong theorem for `S`. `notes/REVISION_LIST.md` R1 to R4 and R8. |
| E10 | the citation of `10.3150/19-BEJ1129` (lines 403, 1372) | The paper describes Theorem 2.3 of that reference as heteroscedastic; the reference (Ding 2020) assumes white noise (Assumption 1.1) and has no variance-profile model. The heteroscedastic input is proved in the tree: `heteroLaw_of_gaussian` (`RMT/Het/Sup.lean`), section 4. Verified 2026-09-02 from arXiv:1702.06975. |
| E11 | `assum:gen_rank_stacksvd_eig_sep` (line 2104), `prop:gen_rank_stacksvd_singleweight` (line 2112) | At a tied root the eigenvectors `z_ℓ` of the tied plane are defined up to a rotation, and the performance formula (a sum of reciprocals) changes with the choice unless the `z_ℓ` are orthogonal for the form `K(γ_ℓ)`; so the distinctness clause matters for the value, not only for the proof. A remark, not an error. Closed form plus Monte Carlo (`scripts/numeric/check_singleweight.py`, seed 20260905). Anchor: `EigSep.sorted` (`RankR/SingleWeight/Scalars.lean`, section 6.8), the field that turns the paper's distinctness clause into a strict order. |
| E12 | the proof of `prop:singleweight_suboptimality` (line 2187) | the branch where the eigenvalue separation assumption fails says only that the performance is bounded above by `β_0²`, with no argument; the Lean proves the exact limit there, the summand of the larger weight in the paper's own display, which is strictly below `β_0²`. A remark, not an error (F18b, 2026-09-07). Anchors: `swLimitEx` (`RankR/SingleWeight/Example.lean`), `hlaw_one` and `prop_singleweight_suboptimality_gaussian` (`RankR/SingleWeight/Suboptimality.lean`). |

Read `notes/paper_edits.md` for the proposed text of each item. Status of every item on 2026-09-02
is `proposed`; none is applied to the paper.

### 7.3 Open questions that concern scope (`notes/FLAGGED.md`)

| # | Question |
|---|---|
| Q2 | When `A_β = I` (exactly one `β_i > 0`) the Python and R reference code disagree and the paper leaves the case open. Lean excludes it through `hthr`. |
| Q3 | Layer 1 is Gaussian-only through `JointGaussianNoise`; `assum:general_noise` is wider. Kept, by user decision 2026-08-29. Since 2026-09-02 (L4) the rank-one SVDstack family takes `IndepNoise` instead; since 2026-09-08 (F31) `thm_theta_est` takes `ThetaEstLaw` (section 2.6), so no rank-one Layer 1 theorem keeps `JointGaussianNoise` (D35, superseded). Since 2026-09-08 (Stage 0 of `notes/NONGAUSSIAN_SCOPE.md`) `ThetaEstLaw` is also discharged for any fixed law with mean 0, variance 1 and finite fourth moment (`thetaEstLaw_of_general`, `thm_theta_est_general`, section 5.4), the first non-Gaussian discharge; the single-table law stays Gaussian-only. |
| Q5 | The `M = 3` remark example sits on the threshold. Now `notes/paper_edits.md` E6. |
| Q6 | `thm:simple_thm1` is cor. 1, not cor. 2. Fixed in `notes/README.md`. |
| Q7 | Five rank-`r` design decisions. Answered: `r_i = 1` first; performance as `tr(Bᵀ specInvTop A r B)`; `SubspaceLaw` as the single new structure; build Ky Fan; the Gaussian rank-`r` law was later proved, so item (5) is closed. |
| Q8 | The paper's removal remark (line 901) is wrong in general. Not formalized. |
| Q9 | The paper writes `prop:stacksvd_general` where it means `prop:stacksvd_subspace` (line 842). Label typo. |
| Q10 | Four appendix results had no Lean: `prop:singleweight_suboptimality`, `prop:gen_rank_stacksvd_singleweight`, `thm:rank_r_stacksvd`, `thm:rank_r_svdstack`. The last two are now proved (sections 6.6, 6.7). The first two now have Lean as well (2026-09-05): both are formalized at Layer 1 (sections 6.8, 6.9); the Gaussian discharge of `SingleWeightLaw` landed the same day (Track G), so `prop:gen_rank_stacksvd_singleweight` has its `_gaussian` theorem, and `prop:singleweight_suboptimality` got its unconditional `_gaussian` theorem on 2026-09-07 (F18b). Only the removal remark remains in the table of section 7.1. |
| Q11, Q12 | The exact heteroscedastic edge. Closed: `heteroEdge_of_gaussian` proves it for Gaussian noise by a sharp Sudakov-Fernique bound, so `thm_stacksvd_weighted_gaussian` is unconditional. |
| Q14 | The second external review of the audit pack. Verdicts recorded. |
| Q15 | The MLE and remark statements. Answered "OK both" on 2026-09-01. |
| Q16 item 1 | Keep the ordered class and state the rank-`r` weighted result at `ℓ_j = j`, or add a permutation `R_i` route for the paper's example. **Open.** Claude's recommendation is to keep the class, since the general-`R_i` theorems already cover the permuted case for svdstack. |
| Q16 item 2 | The component-clause statement form of `thm:rank_r_svdstack` (canonical frame, `hS : StrictAnti`). **Closed**: user OK 2026-09-05 ("prove simplest rank r (S_I sorted)"). |
| Q16 last | Whether `notes/paper_edits.md` E1 to E11 go to the paper before or after the external audit. **Closed** 2026-09-05: after; the audit reads the paper as it is, with E1 to E11 as the list of clarifications the authors will make. |

### 7.4 Numerical checks

Every closed form in the tree was recomputed against an independent implementation. The
reference code is the Python `simulations.py` and the R package `stackedSVD`
(`theory_pred.R`, `getTheoryPred`, `compute_asymptotic_power`, `compute_x_star`) in
the authors' analysis code (the `stackedSVD` repository).

- `notes/archive/audit_numeric_2026-08-29.md` (seed 20260829, printed): `betaSq` against
  `prop:single_table`, the Python `beta_vec` and the R `beta.true^2` (worst 2.2e-16);
  `stackSVDLimit` against `betaSq(‖θ‖₂, ‖c‖₁)`, the Python and the R stack values (worst
  2.2e-16); `svdstackLimit` against the Python and the R svdstack values away from the
  degenerate case (worst 3.3e-16); `rhoSq` against `θ²+1+c+c/θ²` above the threshold and
  `bulkEdge` at or below it (0.0); `cor.2` against `stackSVDLimit` on the selected set
  (2.2e-16). `Abeta` is checked through `svdstackPerfClosed` (2.2e-16), not on its own. The
  one disagreement is the degenerate `A_β = I` case, which Lean excludes through `hthr`
  (question Q2).
- `notes/archive/plan_heterolaw_A.md` and `notes/archive/agent_reports/expert_heterolaw_A.md`:
  the overlap identity `1/(ρ F'(ρ)) = L(w)` at `θ = (1.3, 0.9, 1.6)`, difference `1.6e-11`,
  and `Assumption4 ↔ ρ > bHet` on 1642 draws.
- `notes/archive/audit_het_skeleton_2026-08-31.md`: the heteroscedastic scalars; the same
  overlap identity at three parameter sets to `3.3e-16` or better; the `M = 1`, `w = 1`
  reduction `bHet = bulkEdge c` and `rhoHet = rhoSq`; and the margin coverage figure.
- `notes/archive/audit_independent_TrackE_wave2b_2026-09-02.md`: the concrete-model instantiation of
  `thm_rank_r_stacksvd_gaussian` in Lean (section 6.7 above), and `gammaR` against
  `theory_pred.R` to `1.0e-13`.
- Per-audit scripts under `scripts/`, each with a printed seed and assertions instead of
  eyeballed output. Reruns on 2026-09-02 (numpy only, no Lean):
  `check_remark_examples.py` (seconds, `RESULT: PASS`), `check_rayleigh_bound.py` (seconds,
  `RESULT: PASS`, worst gap `-2.6e-04` against a bound that must be `≤ 1e-12`),
  `check_port_stacksvd_subspace.py` (about 3 minutes, `OVERALL: PASS` on six cells including
  two tie cells and one subcritical cell), `check_edge_count.py` (a Monte Carlo of task U7,
  seed 2026090101; longer than 2 minutes, not rerun in the docs audit).

## 8. Build and verification

```sh
cd lean && lake exe cache get && lake build   # rechecks every proof
cd .. && scripts/check_sorries.sh                   # then the eight gates, see scripts/README.md
```

(On the development server a capped wrapper, `lake-build-capped.sh 12`, stands in for bare `lake build`.)

Build of record: root `lake build` of 2026-09-08 08:14:24 to 08:16:49 EDT
(`lake-build-capped.sh 6`) on the Lean tree of commit `8736b2a` (F31, `thm_theta_est` at
Layer 1), which is that of the current head: **8930 jobs, exit 0**, 0 errors, 3 jobs
rebuilt (`ThetaEst.lean` 39 s, `Main.lean` 34 s, the root 32 s), 5 replayed (the vendored
COLT83 files, the only jobs with a stored log); 0 linter warnings in our 162 files. Seven
gates on that tree, each exit 0 (gates 1 to 6 at 08:30:09 to 08:31:57 EDT, gate 7 at
08:20:53 to 08:29:05 EDT): `check_sorries.sh` 0 sorry sites in 0 declarations; `check_axioms.sh`
**4622 declarations** audited (4622 2329), no axiom outside `propext`, `Classical.choice`
and `Quot.sound`, no untracked `sorry`; `check_root_imports.py` 182 modules, 162 ours, all
162 inside the audited closure; `check_theorems_sigs.py` 128 signatures, 116 verbatim, 12
documented abbreviations; `check_layering.sh` 3 endpoints, 0 forbidden hits;
`check_paper_edits.py` 29 replacements, 6 amended blocks; `check_kernel.sh` 183 modules
replayed through the kernel, 0 problems, every one of the 183 oleans of the package (492
s, 3 workers). Before that (`18452a2`, `Main.lean` and F30, 2026-09-07 21:55:38 to
22:16:22 EDT, 3 rebuilt, 5 replayed; seven gates 22:16:52 to 22:26:48 EDT): 8930 jobs,
4618 declarations (4618 2324), 182 modules, 162 ours, 126 signatures, 114 verbatim, 12
documented abbreviations, 29 replacements, 183 modules replayed (507 s). Before that
(`926701c`, the linter campaign, 2026-09-07 18:14:06 to 18:34:34 EDT, 114 rebuilt, 5
replayed): 8929 jobs, 4576 declarations (4576 2323), 181 modules, 161 ours, 126
signatures, 114 verbatim, 12 documented abbreviations, 29 replacements. Before that (`85e6a0e`, F18c,
2026-09-05 19:23:32 to 19:24:29 EDT, 1 rebuilt, 73 replayed; the header run of Appendix A
of `AUDIT_DOC.md` at 19:36 EDT on the docs commit `9c82538`, the same Lean tree): 8922
jobs, 4519 declarations (4519 2249), 174 modules, 154 ours, 125 signatures, 113 verbatim,
12 documented abbreviations, 28 replacements.

Two gates:

```sh
scripts/check_sorries.sh                        # the tree matches the (empty) docs/SORRIES.md ledger
CHECK_AXIOMS_BUILD=0 scripts/check_axioms.sh    # every declaration: only the 3 standard axioms
```

**Reading `check_axioms.sh`.** The script runs the driver `scripts/AxiomAudit.lean` against
the built oleans. The driver prints one `DECL <name> | <axioms>` line per declaration of the
`StackedSVD` namespace, a `SORRYDEPS` line after any `DECL` that carries `sorryAx`, and a
final `AXIOMAUDIT_END <reported> <skipped>` line. The gate then reads that output and prints
only a count and the failures, so a passing run ends with two lines:

```
check_axioms: 4226 declarations audited (<reported> <skipped>).
check_axioms: OK - no axiom outside {propext, Classical.choice, Quot.sound} and no untracked sorry.
```

Pass `--raw` to see the driver output itself, or `--verbose` to keep the gate and print every
declaration with its axioms. The gate fails when a declaration depends on anything but
`propext`, `Classical.choice` and `Quot.sound`. `sorryAx` is the one conditional case: it
passes only when the declaration has a row in `docs/SORRIES.md`, or inherits the `sorry` from such
a row, and the driver names the source on the `SORRYDEPS` line. `docs/SORRIES.md` is empty, so any
`sorryAx` fails. Without `CHECK_AXIOMS_BUILD=1` the gate reads the oleans, not the sources; it
then prints a freshness report that names every source file with no olean or an olean older
than itself, so a declaration edited since the last build is audited in its old form. Set
`CHECK_AXIOMS_BUILD=1` to build first (it refuses while another `lake` process runs).
`AXIOM_AUDIT_VENDOR=1` also audits `StackedSVD/Vendor/`, which is skipped by default.

`scripts/check_sorries.sh` prints one line, and on this tree it reads

```
check_sorries: OK - 0 sorry sites in 0 declarations, 0 rows, all matched.
```

(run 2026-09-03, exit 0). `--sites` lists every `sorry` with its declaration and `--tracked`
lists every `docs/SORRIES.md` row; both are empty here.

A five-line probe an auditor can paste into a Lean file inside the project:

```lean
import StackedSVD
#print axioms StackedSVD.SpikedModel.singleTableLaw_of_gaussian
#print axioms StackedSVD.MultiTableModel.thm_stacksvd_weighted_gaussian
#print axioms StackedSVD.MultiTableModel.thm_svdstack_weighted_gaussian
#print axioms StackedSVD.UnalignedModelR.thm_rank_r_stacksvd_gaussian
```

Each prints `[propext, Classical.choice, Quot.sound]`.

`scripts/AuditSignatures.lean` prints the elaborated signature and the axioms of every
headline theorem, which is the fastest way to read the statements as Lean sees them.

**A third gate on this document.** `python3 scripts/check_theorems_sigs.py docs/THEOREMS.md`
compares every Lean signature quoted here with the tree and prints a `NO VERBATIM MATCH` line
for each block that does not match character for character, then the totals. That line is
expected for the signatures this document abbreviates on purpose, which section 1 and section
6.4 name; `NAME NOT IN TREE` is always a defect.

**A fourth gate, on coverage.** `python3 scripts/check_root_imports.py` fails when a module of
the package is outside the import closure of `lean/StackedSVD.lean`. It matters because
`AxiomAudit.lean` audits what `import StackedSVD` reaches, so a module missing from that file
would be compiled by `lake build` and never audited. One vendored module,
`Vendor/COLT83/Axioms.lean`, is outside the closure on purpose: it is the vendor's own
`#print axioms` driver that nothing imports.

`scripts/README.md` documents these and every other script, and
`scripts/make_audit_pack_header.sh` runs the build and four of the seven gates in one pass
(`AUDIT_PACK_JOBS` caps the cores of the build it starts; the default is 12).

### 8.1 The hypothesis sets are satisfiable

A theorem whose hypotheses no model meets is true and empty. `lean/StackedSVD/Sat.lean`
settles that question for the three headline Gaussian endpoints. It builds concrete Gaussian
models and applies each endpoint to one of them with no hypothesis left over:

| Theorem | Endpoint it applies | Model |
|---|---|---|
| `Sat.Rank1.sat_stacksvd_weighted` | `thm_stacksvd_weighted_gaussian` | `M = 2`, `c = (1,1)`, `θ = (2,2)` |
| `Sat.Rank1.sat_svdstack_weighted` | `thm_svdstack_weighted_gaussian` | the same model |
| `Sat.RankR.sat_rank_r_stacksvd` | `thm_rank_r_stacksvd_gaussian` | `M = 2`, `r = 2`, `θ = [[3, 3/2], [5/2, 1]]`, `R_i = 1` |

The rank-one construction is parameterized: `n_i N = k_i (N + 1)` and `d N = l (N + 1)`, so the
aspect ratio is exactly `c_i = k_i / l` at every `N` and the regime hypothesis holds by
`div_self`. The noise is `Measure.pi` of `gaussianMatrix`, one factor per table, which gives
`JointGaussianNoise` by `HasLaw.id`. The conclusions are not trivial at these instances: every
model is above the detection threshold, so every limit is strictly positive (rank `r`:
`∑_i θ_i0⁴ / c_i = 120.06` and `∑_i θ_i1⁴ / c_i = 6.06`).

`SingleWeightLaw` (section 6.8) has no separate witness in this module: the natural
candidate is a rank-one model with one table at `θ = 0`, and no
`MultiTableModel`-to-`UnalignedModelR` bridge exists to build it (`SpikedModelR.hθnn` allows
`θ = 0` since F8; the multi-table map is still to write). Since 2026-09-05 the Gaussian
discharge `singleWeightLaw_of_gaussian` (`RankR/SingleWeight/Het/Sup.lean`) produces the
structure on every Gaussian `UnalignedModelR` with `r ≤ ∑ i with w i ≠ 0, n i N` at every
`N` that meets `EigSep`, which shows the hypothesis set consistent whenever such a model with an `EigSep` instance
exists; the two-table witness `mdl` of section 6.9 at `w = (1, 9/10)` is one (the roots
`γ_l = w_l²(1 + θ_0²)` are then `3.56` and `2.88`, both above `wSqMax w = 1`, distinct, with
threshold sums `0.24` and `0.43` below `1` and eigenvectors `e_l`; checked numerically, not
written in Lean).

The module is in the import closure of `StackedSVD.lean`, so `check_axioms.sh` audits the
three witnesses on every run and they cannot rot. The models come from two audits that built
them outside the tree (`notes/archive/audit_independent_D33_2026-09-02.md` and section 5.3 of
`notes/archive/audit_independent_TrackE_wave2b_2026-09-02.md`).

## 9. Record of the cleanup

The branch `cleanup` (2026-09-02) removes duplicate lemmas, drops redundant hypotheses, and
fixes wrong paper labels in docstrings. It adds no result. Every deleted declaration and every
changed signature is listed in `notes/archive/cleanup_2026-09-02_removed.md`. Read that file
before you compare a signature here with an earlier version of the tree.
