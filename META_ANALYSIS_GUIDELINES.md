# Meta-analysis guidelines

> Source: <https://github.com/alexmascension/agent-coding-guidelines> · Version 1.1.1 (2026-10-14)
>
> Key words **MUST**, **SHOULD**, **MAY** follow RFC 2119. Bracketed keys such as `[Borenstein 2023]` point to
> `OMICS_REFERENCES.md`; keys marked with a dagger point to `STATISTICS_REFERENCES.md`.

## 0. How to use this document

* Audience: coding agents and people with introductory statistics. Each rule is followed by its rationale (*Why*, *Trade-off*, *Deviate when*) where it applies.
* Copy this file with `OMICS_REFERENCES.md`, `STATISTICS_GUIDELINES.md` and `STATISTICS_REFERENCES.md`; for omics
  meta-analyses also copy the guideline of the data type (`RNASEQ_BULK_GUIDELINES.md`, `SCRNASEQ_GUIDELINES.md`).
* Part A (§1–8) covers the general methodology of combining results across studies. Part B (§9–13) applies it to
  transcriptomic data, including single-cell data.
* Background texts: the *Cochrane Handbook* and *Introduction to Meta-Analysis* (Borenstein et al.), listed in
  `OMICS_REFERENCES.md` (Books and online resources).

## Part A. General methodology

## 1. Protocol before data

* **Write the protocol before collecting studies**: question, population or tissue, exposure or condition,
  comparison, outcome, inclusion and exclusion criteria, sources searched, effect measure and planned analyses
  (including subgroups and sensitivity analyses). Record it in `PROJECT_SCOPE.md` or the decisions log, with its date.
  *Why:* choosing studies or analyses after seeing results biases the combined estimate
  (`STATISTICS_GUIDELINES.md` §1, design).
  * *Trade-off:* the protocol may need amendments once studies are seen; amendments are allowed if dated and
    justified, and analyses added later are labelled as post hoc.
* **Report the selection process with a PRISMA (Preferred Reporting Items for Systematic reviews and Meta-Analyses) 2020 flow diagram**: records found, screened, excluded with reasons,
  and included [Page 2021].
* Screening and data extraction SHOULD be checked by a second person or a second independent pass, and disagreements
  recorded.

## 2. Data extraction and study quality

* Extract, for every study: the effect estimate, its standard error (or the data to compute it), sample size per
  group, and study-level covariates (design, platform, tissue, population, sex, age).
* **Harmonize the direction and scale of effects** before combining (same reference group, same units, log scale for
  ratios).
  * *Why:* combining a study that reports treated versus control with one that reports control versus treated
    cancels real effects; mixing natural and log scales distorts weights.
* Assess risk of bias per study with a structured tool appropriate to the design (RoB 2, the Cochrane risk-of-bias tool, for randomized trials
  [Sterne 2019]); for omics studies, record at least replication, confounding of batch with condition, and sample
  quality.

## 3. Effect measures

* Use an effect measure that is comparable across studies:
  * continuous outcomes in different units: standardized mean difference (Hedges' $g$);
  * same units: mean difference;
  * binary outcomes: log odds ratio or log risk ratio;
  * expression: $\log_2$ fold change with its standard error.
* Ratios are combined on the log scale and back-transformed only for display.

## 4. Statistical model

* **Default: random-effects model.** Studies differ in population, protocol and platform, so their true effects are
  expected to differ. The random-effects mean is the average of a distribution of true effects
  [DerSimonian 1986; Riley 2011].
  * *Trade-off:* random-effects weights are more equal across studies than fixed-effect weights, so small studies
    gain influence; small-study effects (§6) then matter more.
* A fixed-effect (common-effect) model is justified only when the studies estimate the same true effect (same
  population and protocol); state the justification.
* **Estimate the between-study variance $\tau^2$ with REML (restricted maximum likelihood) or Paule-Mandel, not DerSimonian-Laird**
  [Veroniki 2015; Langan 2018].
  * *Why:* DerSimonian-Laird underestimates $\tau^2$ when heterogeneity is substantial or studies are few, which
    makes confidence intervals too narrow [Veroniki 2015; Langan 2018].
* **Compute the confidence interval of the mean with the Hartung-Knapp-Sidik-Jonkman (HKSJ) method**
  [IntHout 2014; Langan 2018]. With very few studies, HKSJ can be anticonservative or erratic; use the modified
  version [Röver 2015] or report both.
  * *Why:* the standard (Wald) interval ignores the uncertainty in $\tau^2$ and is too narrow with few or unequally
    sized studies; HKSJ uses a $t$ distribution with an adjusted variance [IntHout 2014].
* **With fewer than about 5 studies, $\tau^2$ is poorly estimated.** Report this limitation; a Bayesian
  random-effects model with an informative prior on $\tau$ MAY be used.

## 5. Heterogeneity

* Report together: $\tau^2$ (or $\tau$), $I^2$, Cochran's $Q$ with its p-value, **and the prediction interval**
  [Higgins 2002; Higgins 2003; Borenstein 2023].
  *Why:* $I^2$ is the proportion of observed variance due to true differences, not the size of those differences; the
  prediction interval gives the range in which the true effect of a new study is expected [Riley 2011].
* The 95% prediction interval of a random-effects meta-analysis with $k$ studies is
  $$\hat\mu \pm t_{k-2,\,0.975}\sqrt{\hat\tau^2 + \widehat{SE}(\hat\mu)^2}$$
* Explore heterogeneity with subgroup analysis or meta-regression only for covariates planned in the protocol
  (§1, protocol), and only with enough studies (rule of thumb: about 10 per covariate).

## 6. Small-study effects and publication bias

* Draw a funnel plot and run asymmetry tests (Egger) only with at least 10 studies [Sterne 2011; Egger 1997].
  * *Why:* with fewer studies the tests have too little power to distinguish asymmetry from chance [Sterne 2011].
* **Asymmetry is not proof of publication bias**: heterogeneity, chance and small-study effects also produce it. Report
  "risk of publication bias" rather than its presence [Afonso 2023].
* Prefer sensitivity analyses that quantify how strong selective publication would have to be to change the
  conclusion [Mathur 2020], or model-averaged adjustment (robust Bayesian meta-analysis) [Bartoš 2022].
* When several effect sizes come from the same study, bias tests and models MUST account for this dependence
  (multilevel models) [Nakagawa 2021].

## 7. Sensitivity analyses

* Leave-one-out analysis and influence diagnostics.
* Exclusion of studies at high risk of bias.
* Alternative $\tau^2$ estimators and confidence-interval methods (§4, statistical model).

## 8. Reporting

* Forest plot with study estimates, confidence intervals, weights, the combined estimate, its confidence interval and
  the prediction interval.
* A table of all included studies with extracted data and covariates, and the PRISMA checklist [Page 2021].
* Software and version; the R package metafor is the reference implementation [Viechtbauer 2010; Balduzzi 2019].

## Part B. Transcriptomic meta-analysis

## 9. Choosing the strategy

Decide, in the protocol, which of these strategies is used and why [Ramasamy 2008; Tseng 2012; Toro-Domínguez 2020]:

| Strategy | When | Main risk |
|---|---|---|
| **Pooled re-analysis**: uniformly reprocess raw data and fit one model with `study` as covariate | Raw data available, same technology | Study effects confounded with condition when a study contains only one group |
| **Effect-size combination**: per-study differential expression (DE), then random-effects meta-analysis of $\log_2$FC and SE per gene | Default for heterogeneous studies | Needs comparable effect measures and standard errors |
| **p-value combination** (Fisher, Stouffer) | Effect sizes not comparable across platforms | Loses direction and magnitude unless done per direction |
| **Rank-based** (RankProd, SumRank) | Mixed platforms, outlier-prone studies | Less interpretable effect sizes |

* **Effect-size combination is the default.** Reproducibility improves more with more datasets than with more samples
  per dataset, and stringent effect-size thresholds with relaxed significance cut-offs gave more reproducible genes
  [Sweeney 2016].
  * *Why:* combining effect sizes keeps direction and magnitude and allows heterogeneity to be estimated per gene;
    p-value combination only says that some study shows an effect.
* p-value combination MUST be done separately for each direction (one-sided p-values), so that up- and
  down-regulation in different studies do not combine into a significant result [Rau 2014]. For dependent tests, the
  harmonic mean p-value controls the family-wise error rate [Wilson 2019].
* Rank products are an option for microarray and mixed-platform data [Hong 2006].

## 10. Data sources and harmonization

* **Reprocess raw data with one pipeline, genome build and annotation release** when possible
  (`RNASEQ_BULK_GUIDELINES.md` §2–4, reference and quantification). Uniformly processed resources: recount3
  [Wilks 2021], ARCHS4 [Lachmann 2018] for human and mouse RNA-seq; CZ CELLxGENE Discover and Census for single-cell
  data [CZI 2024].
  * *Why:* differences in reference, annotation and quantifier produce study-specific gene-level differences
    (`RNASEQ_BULK_GUIDELINES.md` §2) that a meta-analysis would interpret as heterogeneity. *Trade-off:*
    reprocessing costs compute; when raw data are unavailable, record each study's pipeline as a study covariate.
* Map genes to a common identifier set and report genes missing from each study; a gene absent from a study is
  missing, not zero.
* Curate study-level metadata into one table with controlled vocabulary (`CODING_GUIDELINES.md` §6).
* Prefer full per-gene summary statistics over published lists of significant genes: lists are already filtered by
  significance and bias the combination.

## 11. Single-cell meta-analysis

* **Per study, run pseudobulk DE per cell type** (`SCRNASEQ_GUIDELINES.md` §10, condition DE), then combine across
  studies. MUST NOT merge cells from several studies and run cell-level tests.
  * *Why:* merging cells across studies confounds study with condition and reintroduces pseudoreplication at a
    larger scale; per-study pseudobulk keeps the donor as the unit and lets heterogeneity between studies be
    measured.
* **Harmonize cell-type labels** across studies with a common reference (label transfer or CellTypist), and check that
  each cell type is replicable across studies with MetaNeighbor [Crow 2018].
* DE genes from single studies of neurodegenerative diseases reproduce poorly across studies; a non-parametric
  reproducibility-based method (SumRank) identified DE genes with higher specificity and sensitivity than merging
  datasets or aggregating p-values [Nakatsuka 2025].
* Batch-corrected values MUST NOT be used for DE across studies; modelling study or batch as a covariate on
  uncorrected data performs better [Nguyen 2023].

## 12. Reporting (omics)

1. Table of datasets: accession, technology, tissue, species, $n$ per group, processing, inclusion reason.
2. Per gene: combined estimate, confidence interval, $\tau^2$, $I^2$, number of studies contributing, and direction
   consistency (number of studies with the same sign).
3. Multiple-testing correction across genes of the combined p-values (`STATISTICS_GUIDELINES.md` §7, multiple
   testing).
4. Sensitivity of the main results to leaving each study out.

## 13. Python implementation

| Task | Python | R reference (via script) |
|---|---|---|
| Random-effects meta-analysis (Paule-Mandel, DerSimonian-Laird) | `statsmodels.stats.meta_analysis.combine_effects` | `metafor::rma` |
| HKSJ confidence interval | `combine_effects(..., use_t=True)`; row `"random effect wls"` of `summary_frame()` | `metafor::rma(test = "knha")` |
| REML $\tau^2$, meta-regression | `pymare` | `metafor::rma(method = "REML")` |
| Prediction interval | compute with the formula of §5 (heterogeneity) | `metafor::predict` |
| p-value combination | `scipy.stats.combine_pvalues` (Fisher, Stouffer) | `metaRNASeq` |
| Publication-bias sensitivity | — | `PublicationBias`, `RoBMA` |
| Cell-type replicability | — | `MetaNeighbor` |
| Public data | — | `recount3`; ARCHS4 and CELLxGENE Census via their Python clients |
