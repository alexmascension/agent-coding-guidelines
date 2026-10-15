# Cell-type deconvolution guidelines

> Source: <https://github.com/alexmascension/agent-coding-guidelines> · Version 2.0.0 (2026-10-15)
>
> Key words **MUST**, **SHOULD**, **MAY** follow RFC 2119. Bracketed keys point to `OMICS_REFERENCES.md`; keys with a
> dagger point to `STATISTICS_REFERENCES.md`. Each rule is followed by ***Why***, ***Trade-off*** and
> ***Deviate when*** where they apply.
> A rule without a key word is a default with the force of SHOULD. Rules carry stable identifiers (e.g.
> `DECON-EST-001`); overrides refer to them (`CODING_GUIDELINES.md` §0, rule identifiers).

## 0. Scope

* Estimating cell-type composition of bulk RNA-seq samples from a single-cell or single-nucleus reference, and using
  the estimates downstream. Spot deconvolution of spatial data is covered in `SPATIAL_GUIDELINES.md` §7; §8 below lists
  what transfers.
* Copy with `RNASEQ_BULK_GUIDELINES.md`, `SCRNASEQ_GUIDELINES.md` and `OMICS_REFERENCES.md`.
* Main use case in mind: postmortem human brain and spinal cord bulk RNA-seq with snRNA-seq references; the rules are
  general.

## 1. What deconvolution estimates

* `DECON-EST-001` **The estimand depends on the method; record it.** Without a correction for RNA content per cell type, an
  estimate is closer to the fraction of RNA contributed by each cell type than to the fraction of cells. Some
  methods include a correction (MuSiC uses cell-size factors estimated from the reference [Wang 2019]); others
  estimate fractions of reads or report cell proportions relative to the reference cells. RNA fractions, cell
  proportions and abundance scores are not synonyms.
  * *Why:* a neuron contributes far more RNA than a microglial cell. In human brain, common methods gave concordant but
    biased neuronal fractions, and even a matched snRNA-seq reference did not remove the bias [Sosina 2021]. Absolute
    cell fractions require mRNA-content renormalization [Racle 2017; Monaco 2019].
  * *Consequence:* name the estimand in results ("estimated RNA fraction", "estimated cell proportion"), treat
    absolute values as uncertain unless validated against cell counts (§5), and prefer comparisons between groups.
* `DECON-EST-002` **Method contract table**: for each method used, the project records the input scale of bulk and reference (as
  documented by the method), the genes used, whether cell size or RNA content is modelled, the estimand, the
  constraints (non-negative, sum to 1) and the validation result. Examples, to be checked against the documentation of
  the installed version:

  | Method | Input scale (documented) | RNA content / cell size | Reported quantity |
  |---|---|---|---|
  | Bisque | bulk read counts; single-cell counts | not modelled explicitly | cell-type proportions relative to the reference |
  | MuSiC / MuSiC2 | linear scale | cell-size factors from the reference | cell-type proportions |
  | BayesPrism / InstaPrism | raw counts (bulk and reference) | not modelled | fraction of reads per cell type; cell-type expression |
  | hspe | log-scale expression | not modelled | mixing proportions |
  | NNLS baseline (§4) | linear, same normalization as the reference profiles | not modelled | RNA fractions |

* `DECON-EST-003` **Accuracy depends on how distinct the cell types are**: methods work for cell types with specific signatures and
  fail for "fuzzy" ones [Sturm 2019]. Closely related subtypes are estimated less reliably than major classes.

## 2. Reference

* `DECON-REFR-001` **Same tissue and, ideally, same region** as the bulk samples. Brain signatures are region-specific, and signatures
  from cultured cells perform worse [Sutton 2022].
* `DECON-REFR-002` **Contains all cell types expected in the tissue**; missing types degrade every estimate [Avila Cobos 2020].
* `DECON-REFR-003` **Several donors**: methods such as MuSiC weight genes by their consistency across subjects [Wang 2019], and
  multi-reference ensembles reduce study-specific bias [Dong 2020].
* `DECON-REFR-004` **Granularity matches the question**: start with major classes (neurons by broad class, astrocytes,
  oligodendrocytes, OPCs, microglia, endothelial cells); go finer only if validation (§5) supports it.
* `DECON-REFR-005` **Nuclei versus whole-tissue bulk**: snRNA-seq lacks cytoplasmic RNA, while tissue bulk contains it. Prefer methods
  that model technical differences between reference and bulk (Bisque) [Jew 2020], or a reference profiled on the same
  material [Avila Cobos 2023].
* `DECON-REFR-006` **Marker selection**: select genes specific to each cell type in the reference (e.g. ratio of expression in the
  target type to the next highest type, the "Mean Ratio" approach) [Huuki-Myers 2025].
* `DECON-REFR-007` **Reference QC (quality control)**: ambient-corrected counts, doublets removed, annotation verified (`SCRNASEQ_GUIDELINES.md` §4 and §9).

## 3. Data preparation

* `DECON-PREP-001` **Use the input scale that each method documents.** Linear mixing models (NNLS, MuSiC, Bisque, BayesPrism) express
  a bulk profile as a weighted sum of cell-type profiles, which is additive only on the linear scale
  [Avila Cobos 2020]; some methods (hspe) are designed for log-scale input.
* `DECON-PREP-002` Use the same gene identifiers and annotation release in reference and bulk; keep genes measured in both.
* `DECON-PREP-003` Normalize only as the method requires: some take raw counts and normalize internally (Bisque, BayesPrism); for
  methods without their own normalization, use counts per million (CPM) or transcripts per million (TPM)
  consistently in bulk and reference. Never use VST (variance-stabilizing transformation) values, and never
  log-transformed values for a method that expects linear input.

## 4. Methods

* `DECON-METH-001` **Default for brain tissue with an snRNA-seq reference: Bisque** (R, via script).
  * *Why:* in a postmortem prefrontal cortex benchmark validated against RNAScope/immunofluorescence, Bisque and hspe
    were the most accurate [Huuki-Myers 2025]; Bisque is designed for technical differences between single-cell
    references and bulk [Jew 2020].
* `DECON-METH-002` **Run a second, methodologically different method** and a simple baseline:
  * MuSiC / MuSiC2 (multi-subject; MuSiC2 when bulk conditions differ from the reference) [Wang 2019; Fan 2022];
  * BayesPrism, which also estimates cell-type-specific expression, or its fast implementation InstaPrism
    [Chu 2022; Hu 2024];
  * non-negative least squares (NNLS) on marker genes as a baseline.
  * *Why:* disagreement across methods with different assumptions shows that the estimate depends on the model.
    Agreement is necessary but not sufficient: methods can share a bias (e.g. RNA content, an incomplete reference). In spatial benchmarks a simple regression baseline outperformed almost half of dedicated methods
    [Sang-aram 2024], so a method that does not beat NNLS adds nothing.
* `DECON-METH-003` Other methods: CIBERSORTx [Newman 2019] (successor of CIBERSORT [Newman 2015]), DWLS [Tsoucas 2019], SCDC [Dong 2020], SQUID [Avila Cobos 2023], deep
  learning (Scaden, TAPE; Python) [Menden 2020; Chen 2022]. Use them when a benchmark on comparable data supports them.
* `DECON-METH-004` Cross-platform settings (reference and bulk from different technologies) change rankings: SpatialDecon and
  cell2location were the most consistent across platforms in one benchmark [Singh 2026].
* `DECON-METH-005` There is no single best method across tissues and conditions [Jin 2021; Sturm 2019]; the choice is justified by
  benchmarks on similar tissue or by local validation.

## 5. Validation

* `DECON-VAL-001` **Pseudobulk validation**: build pseudobulk samples from reference donors **not used** to build the signature, with
  known composition, and compare estimates with the truth: correlation, absolute error and bias (mean signed error)
  per cell type, over the range of abundances expected in the bulk samples (including rare types), and with an
  incomplete reference (one cell type removed) to see where its signal goes.
  * *Why:* using the same cells for signature and test overestimates accuracy (circularity,
    `STATISTICS_GUIDELINES.md` §6).
* `DECON-VAL-002` **Orthogonal measurements** where available: histology, immunostaining, snRNA-seq proportions of the same samples.
* `DECON-VAL-003` **Sanity checks**: proportions are non-negative and sum to 1 (or are explicitly unconstrained); known biology is
  recovered (e.g. higher oligodendrocyte fraction in white-matter-rich samples; neuronal loss in disease where
  expected); estimates are stable when one reference donor is removed.
* `DECON-VAL-004` **Report the uncertainty**: method agreement, validation error per cell type, and cell types that failed
  validation, which MUST NOT be interpreted.

## 6. Using the estimates downstream

* `DECON-DOWN-001` **As covariates in DE (differential expression)**, estimated proportions reduce false associations caused by composition differences
  [Patrick 2020].
  * *Trade-off:* when composition change is the disease effect (e.g. neuronal loss in ALS), adjusting for it changes
    the question from "what differs in the tissue" to "what differs given the composition". Report both analyses and
    state which question each answers.
  * *Collinearity:* proportions sum to 1; include all but one cell type, or transform (e.g. isometric log-ratio), and
    avoid including many highly correlated types.
* `DECON-DOWN-002` **Comparing proportions between groups**: treat them as compositional data (an increase in one type lowers the
  others) and use the donor as the unit (`STATISTICS_GUIDELINES.md` §3.2, data type to model).
* `DECON-DOWN-003` **Cell-type-specific expression** inferred from bulk (BayesPrism, CIBERSORTx) is a model estimate with large
  uncertainty for minor cell types; validate with single-cell data before interpreting.

## 7. What to report

1. `DECON-REP-001` Reference: source, tissue and region, number of donors and cells, cell-type annotation and granularity, processing.
2. `DECON-REP-002` Marker selection method and number of genes per cell type.
3. `DECON-REP-003` Bulk input scale and normalization.
4. `DECON-REP-004` Methods, versions and parameters; baseline results.
5. `DECON-REP-005` Validation design and error per cell type; cell types excluded from interpretation.
6. `DECON-REP-006` Downstream use: whether DE was run with and without proportions, and the transformation used.

## 8. Relation to spatial spot deconvolution

* `DECON-SPAT-001` Transfers directly: complete, matched reference; documented input scale; validation with held-out data; the
  estimand depends on the method (§1).
* `DECON-SPAT-002` Differs: spots contain few cells, so sparsity and spatial smoothing matter; spatial methods (cell2location, RCTD)
  model platform effects and per-spot detection [Kleshchevnikov 2022; Cable 2021]; see `SPATIAL_GUIDELINES.md` §7.

## 9. Implementation

| Method | Language | Notes |
|---|---|---|
| Bisque | R (`BisqueRNA`) | Default for brain with snRNA reference |
| hspe | R | Top performer with Bisque in postmortem cortex [Huuki-Myers 2025] |
| MuSiC / MuSiC2 | R | Multi-subject reference |
| BayesPrism / InstaPrism | R | `pybayesprism` exists in Python but is at an early version (0.1) |
| NNLS baseline | Python (`scipy.optimize.nnls`) | On marker genes, linear scale |
| Scaden | Python (`scaden`) | Deep learning; needs simulated training data |
| Marker selection (Mean Ratio) | R (`DeconvoBuddies`) | Or reimplement the ratio in Python with a test against the R output |
