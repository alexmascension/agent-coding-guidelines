# Statistics guidelines

> Source: <https://github.com/alexmascension/agent-coding-guidelines> · Version 2.0.0 (2026-10-15)
>
> Key words **MUST**, **MUST NOT**, **SHOULD**, **SHOULD NOT** and **MAY** follow RFC 2119, as in
> `CODING_GUIDELINES.md`. Bracketed keys such as `[Squair 2021]` point to entries in
> `STATISTICS_REFERENCES.md`.
> A rule without a key word is a default with the force of SHOULD. Rules carry stable identifiers (e.g.
> `STAT-DES-001`); overrides refer to them (`CODING_GUIDELINES.md` §0, rule identifiers).

## 0. How to use this document

* `STAT-USE-001` Audience: coding agents and people with introductory statistics (one or two university
  courses). Each rule is followed by its rationale so the reasoning can be checked, not only
  followed: *Why* (mechanism), *Trade-off* (cost) and *Deviate when* (justified exceptions).
* `STAT-USE-002` Copy this file and `STATISTICS_REFERENCES.md` verbatim into the project root when the project
  involves statistical analysis. Project-specific deviations go in the *Overrides* section of
  `PROJECT_CONTEXT.md`, as described in `CODING_GUIDELINES.md` §0 (how to use the guidelines).
* `STAT-USE-003` This file covers principles that apply to any data type. Tool-specific procedures for
  sequencing data (normalization, differential expression (DE) packages, single-cell workflows) belong in
  `RNASEQ_BULK_GUIDELINES.md`, `SCRNASEQ_GUIDELINES.md` and `SPATIAL_GUIDELINES.md`; figure design belongs in `FIGURE_GUIDELINES.md`.
* `STAT-USE-004` Priority when rules conflict: **validity of the inference > transparency of what was done >
  statistical power > convenience**.

## 1. Design before analysis

* `STAT-DES-001` **Write down the question as a comparison before touching the data**: what is compared with
  what, in which units, and which result would count as an effect. Record it in
  `PROJECT_SCOPE.md` or the decisions log.
  *Why:* an analysis chosen after seeing the data can always find something; this inflates
  false positives [Simmons 2011].
  * *Deviate when:* exploratory analyses are allowed, but are labelled as exploratory, and their findings are
    treated as hypotheses to confirm in independent data.
* `STAT-DES-002` **Identify the experimental unit**: the smallest entity that was independently assigned to a
  condition (a mouse, a patient, a culture flask). The sample size $n$ is the number of
  experimental units, not the number of measurements [Lazic 2018].
  *Example:* 3 mice per group, 500 cells measured per mouse: $n = 3$ per group, not 1500.
  * *Trade-off:* the honest $n$ is often small, which lowers power; the remedy is more units, not counting
    measurements as units.
* `STAT-DES-003` **Distinguish biological from technical replicates.** Biological replicates (independent
  units) support conclusions about the population; technical replicates (repeated measurement of
  the same unit) only measure assay noise and are averaged or modelled, never counted as $n$
  [Blainey 2014].
* `STAT-DES-004` **Batch MUST NOT be confounded with the condition of interest.** Spread every condition across
  processing days, plates, lanes and operators, and record these variables in the sample
  metadata [Leek 2010].
  *Why:* if all controls were sequenced on Monday and all treated samples on Tuesday, the two
  effects cannot be separated by any statistical method.
  * *Trade-off:* balancing requires planning sample processing in advance. *Scope:* the rule applies where processing can be planned; when balancing is impossible
    (e.g. retrospective samples); then record the confounding and report it as a limitation, because no analysis
    can remove it.
* `STAT-DES-005` **Justify the sample size** and state how: available resources, an a-priori power
  calculation, a target precision for the estimate, or field conventions [Lakens 2022]. For bulk
  RNA-seq, fewer than 6 biological replicates per group leaves many true changes undetected
  [Schurch 2016].
  * *Why:* an unjustified sample size usually hides low power, and with low power the significant results that do
    appear overestimate the effect (last rule of this section).
* `STAT-DES-006` **MUST NOT report post-hoc ("observed") power.** It is a transformation of the p-value and adds
  no information [Hoenig 2001].
* `STAT-DES-007` **Plan positive and negative controls** (a gene or condition with a known effect, and one with
  none) to check that the pipeline detects what it should and nothing else.
* `STAT-DES-008` Low power has two consequences beyond missed effects: for a given significance level and proportion of true
  effects among the hypotheses tested, significant results are more likely to be false, and their effect sizes are
  exaggerated (type M error) or even have the wrong sign
  (type S error) [Button 2013; Gelman & Carlin 2014].

## 2. Exploration and quality control

* `STAT-EXPL-001` **Plot the raw data before modelling**: distributions per sample, scatter plots of key
  variables, PCA (principal component analysis) of samples coloured by condition and by every batch variable.
  *Why:* outlier samples, label swaps and batch structure are usually visible before any test.
* `STAT-EXPL-002` **Missing is not zero.** Decide for every variable what a missing value means (not measured,
  below detection, not applicable) and encode it as null, never silently as 0.
  * *Why:* when values are missing more often in one group (e.g. below detection in diseased samples), encoding
    them as 0 creates a difference that reflects missingness, not biology.
* `STAT-EXPL-003` **Zeros must be explicit.** When counting events per group (e.g. cells of a type per sample),
  combinations with no events MUST appear as rows with count 0; otherwise means and tests are
  computed on the wrong denominator.
* `STAT-EXPL-004` **Outliers are removed only by a rule fixed before looking at the results** (e.g. a quality control (QC) metric
  threshold), applied blindly to all groups, and reported with the number of units removed.
  MUST NOT remove points because they weaken an effect.
  * *Trade-off:* a fixed rule may keep a few bad points or remove a few good ones; that is preferable to a decision
    that depends on the result. *Deviate when:* a measurement is demonstrably erroneous (documented technical
    failure); remove it and record the reason.
* `STAT-EXPL-005` **MUST NOT transform the response or delete data to force normality.** Choose a model that fits
  the data type instead (§3, model choice).
* `STAT-EXPL-006` **Check implementation details on a few cases by hand** (log fold change with its pseudocount,
  normalization direction, default filters). Package documentation and behaviour sometimes
  differ.

## 3. Choosing the model

### 3.1 Principle: tests are models

* `STAT-MOD-001` **Prefer fitting a model and testing on it over standalone tests.** Most classical tests are
  special cases of linear or generalized linear models: the Student (equal-variance) two-sample t-test is an
  ordinary least-squares regression on a binary group variable, while the Welch t-test corresponds to a model with a
  separate variance per group; one-way ANOVA (analysis of variance) is a regression on a categorical variable; the
  Mann-Whitney test is the score test of a proportional-odds model with one binary predictor (the test coincides; the
  estimands and default reports of the two are not interchangeable).
  *Why:* a model can include covariates (age, sex, batch), interactions and arbitrary contrasts,
  and makes the assumptions explicit (Practitioner sources, Olszewski).
  * *Trade-off:* models must be specified correctly and can fail to converge; check convergence and residuals.
    *Deviate when:* a simple two-group comparison without covariates; the Welch t-test is already the right model.
* `STAT-MOD-002` **The assumptions concern the residuals, not the raw data.** For a linear model: residuals
  independent, with roughly constant variance, and (for small $n$) roughly normal. Check them with
  residual plots, not with a normality test of the raw values.
* `STAT-MOD-003` **When an assumption fails, change the model, not the test family.** See the table in §3.2
  (data type to model).

### 3.2 Data type to model

| Response | Default model | Notes |
|---|---|---|
| Continuous, roughly symmetric | Linear model (OLS, ordinary least squares); Welch t-test for two groups | Use heteroscedasticity-robust errors (`HC3`) if variances differ. |
| Continuous, positive and right-skewed | Gamma GLM (generalized linear model; log link) or linear model on log scale | State which; they estimate different quantities (mean vs geometric mean). |
| Counts | Negative binomial GLM; Poisson only if no overdispersion | Sequencing counts: use dedicated tools (`RNASEQ_BULK_GUIDELINES.md`). |
| Binary (yes/no) | Logistic regression | Reports odds ratios; a classifier is a separate thresholding step. |
| Proportion out of $k$ trials | Binomial GLM with $k$ as trials | Do not analyse percentages as continuous values when $k$ is known. |
| Proportion without known $k$, strictly between 0 and 1 | Beta regression | Its support is (0, 1): exact 0s or 1s need a zero/one-inflated beta model or a documented transformation. With known $k$ and extra variation, use a beta-binomial or quasi-binomial model. |
| Ordered categories (scores, grades) | Proportional-odds (ordinal logistic) model | |
| Time to event | Cox proportional hazards or parametric survival model | Censored observations MUST NOT be dropped. |
| Compositions (relative abundances summing to 1) | Log-ratio methods (CLR: centred log-ratio, ALR: additive log-ratio) | Correlations between raw proportions are spurious [Gloor 2017]. |
| Repeated or nested measurements | Mixed model, GEE: generalized estimating equations, or aggregation to the experimental unit | §4 (hierarchical data). |

### 3.3 Specific rules

* `STAT-MOD-004` **Two-group comparison of means: Welch's t-test by default**, not Student's.
  *Why:* Welch does not assume equal variances and loses almost nothing when they are equal
  [Delacre 2017].
  * *Trade-off:* a negligible loss of power when variances are truly equal. *Deviate when:* the data are paired
    (use a paired test or a model with a pairing term).
* `STAT-MOD-005` **MUST NOT choose the test by first testing normality** (e.g. Shapiro-Wilk, then t-test or
  Mann-Whitney depending on the result). Choose the method a priori from the data type and the
  question.
  *Why:* the two-step procedure distorts the error rates of the final test [Rochon 2012]; with
  large samples the t-test is valid for non-normal data anyway [Fagerland 2012; Lumley 2002].
* `STAT-MOD-006` **Rank tests answer a different question.** Mann-Whitney tests whether a value from one group tends to be
  larger than a value from the other, i.e. whether $P(X > Y) + \tfrac{1}{2}P(X = Y)$ differs from 1/2; Kruskal-Wallis
  extends this to several independent groups (at least one group tends to differ), and Friedman to blocked or
  repeated-measures designs. State the hypothesis in those terms, with the dependence structure. They MUST NOT be
  reported as comparisons of medians unless the distributions differ only by a shift [Divine 2018].
  * *Deviate when:* the outcome is ordinal or has extreme outliers and the question is about which group tends to
    have larger values; report the result in those terms (and, under a shift model, a Hodges-Lehmann shift
    estimate).
* `STAT-MOD-007` **MUST NOT dichotomize continuous variables** (e.g. "high vs low expression" at the median)
  for testing. Model the variable as continuous, with splines if the relation is non-linear.
  *Why:* dichotomizing discards information, loses power and can create spurious effects
  [Altman 2006].
  * *Scope:* the rule does not cover questions where a clinically defined cut-off is itself the question (e.g. a diagnostic threshold); state the
    cut-off and its source.
* `STAT-MOD-008` **Randomized two-group before/after designs: adjust for the baseline (ANCOVA, analysis of covariance)**, i.e.
  model the follow-up value with the baseline as covariate, rather than analysing change scores or comparing
  within-group p-values [Vickers 2001].
  * *Scope:* the argument relies on randomization (groups equal at baseline in expectation). In non-randomized
    comparisons, single-group designs or repeated measures over time, the choice of model (mixed model, change
    scores, ANCOVA) depends on the estimand and on how groups were formed; state it.
* `STAT-MOD-009` **Test interactions directly.** "The effect is significant in males and not in females" is not
  evidence that the effect differs by sex; fit the interaction term `treatment:sex` and test it
  [Gelman & Stern 2006; Nieuwenhuis 2011].
  * *Trade-off:* interaction effects are estimated with less precision than main effects, so a study powered for
    main effects is often underpowered for interactions; report the interaction estimate and its interval, not only
    its p-value.
* `STAT-MOD-010` **Wald, likelihood-ratio (LRT) and score tests** usually agree. When they disagree markedly
  (small samples, estimates near a boundary, separation in logistic regression), the large-sample
  approximation is unreliable for all of them, the LRT included. Diagnose the cause and use a method designed for
  it: penalized likelihood (Firth) for separation, parametric bootstrap or boundary-corrected reference
  distributions for variance components, exact methods for small counts; report the disagreement.
* `STAT-MOD-011` **Encode categorical predictors explicitly** and state the reference level; coefficients are
  differences with respect to it.
  * *Why:* software picks the reference level alphabetically by default; a change in level names silently changes
    the meaning of every coefficient.

## 4. Hierarchical data and pseudoreplication

* `STAT-HIER-001` **Pseudoreplication** is treating measurements that are not independent as if they were (cells
  from the same mouse, repeated time points of one patient, technical replicates). It MUST be
  avoided [Lazic 2010].
  *Why:* the test believes it has far more information than it has; standard errors shrink and
  false positives multiply. Single-cell DE methods that treat cells as replicates report hundreds
  of false DE genes [Squair 2021].
* `STAT-HIER-002` Use one of three valid approaches, in order of simplicity:
  1. **Aggregate to the experimental unit** (mean per animal, pseudobulk counts per sample and
     cell type) and analyse the aggregates.
  2. **Mixed model** with a random effect for the unit (`(1 | mouse)`) [Aarts 2014; Zimmerman
     2021].
  3. **GEE** with the unit as cluster, when population-averaged effects are the target.
* `STAT-HIER-003` For confirmatory tests with mixed models, include the random effects implied by the design
  (random slopes when the condition varies within units); too simple a structure inflates type I
  error [Barr 2013].

## 5. Batch effects and technical covariates

* `STAT-BATCH-001` **Batch goes into the model as a covariate** (`~ batch + condition`), not into a pre-corrected
  data matrix that is then tested.
  *Why:* a pre-corrected matrix looks like it has more degrees of freedom than it has, and with
  unbalanced designs correction that "preserves group differences" exaggerates them
  [Nygaard 2015].
  * *Trade-off:* each batch level uses a degree of freedom, which matters with few samples per batch.
* `STAT-BATCH-002` **Batch-corrected matrices are for visualization** (PCA, heatmaps, clustering) or for methods
  that cannot model covariates; state this wherever they are used.
* `STAT-BATCH-003` **Verify from metadata what a "batch" is.** A replicate number is a batch only if it reflects a
  shared technical event (same day, plate, run) across conditions; independent biological
  replicates labelled 1, 2, 3 are not a batch (Practitioner sources, Tang).
* `STAT-BATCH-004` If batch and condition are completely confounded (§1, design), no correction is possible:
  report it as a limitation instead of correcting.

## 6. Circular analysis and data leakage

* `STAT-CIRC-001` **MUST NOT use the same data to select and to test** without a method that accounts for the
  selection [Kriegeskorte 2009]. Common forms:
  * Clustering cells and then testing for markers between those clusters: the p-values are too
    small because the clusters were built to differ [Gao 2022]. Treat them as ranking scores, or
    use count splitting [Neufeld 2022] or a synthetic null [Lee 2025].
  * Selecting the genes most correlated with an outcome and then reporting their correlation.
  * Defining a region of interest from the effect and then testing the effect inside it.
* `STAT-CIRC-002` **Machine learning: every data-dependent step goes inside the cross-validation loop**, fitted
  on training folds only: feature selection, scaling, imputation, batch correction, hyperparameter
  tuning (nested cross-validation) [Ambroise 2002; Varma 2006; Kapoor 2023].
  * *Trade-off:* nested cross-validation multiplies computation; it is the cost of an unbiased error estimate.
* `STAT-CIRC-003` **The test set is used once**, after all choices are fixed.
  * *Why:* every look at the test set followed by a change in the model is a fit to the test set, and the final
    error estimate becomes optimistic.
* `STAT-CIRC-004` **MUST NOT draw quantitative conclusions from 2-D embeddings** (t-SNE: t-distributed stochastic neighbour embedding, UMAP: uniform manifold approximation and projection): distances and
  cluster sizes in the plot are not interpretable [Chari 2023]. Use them for display only.

## 7. Multiple testing

* `STAT-MT-001` **Define the family of tests** before correcting: all genes in one contrast, all contrasts of
  one model, all regions… State the family size $m$ with the results.
* `STAT-MT-002` **Default: Benjamini-Hochberg (BH) adjusted p-values**, controlling the false discovery rate
  (FDR): among the results called significant, the expected fraction of false positives is at
  most the chosen level (e.g. 5%) [Benjamini 1995]. The FDR is an expectation over repetitions of the study; the
  fraction of false positives in one particular list can be higher.
  * *Trade-off:* FDR control accepts that some reported results are false (e.g. up to 5%); it is suited to
    discovery lists, not to claims about single results.
* `STAT-MT-003` **Use family-wise control (Holm; Bonferroni) when any single false positive is costly** (few
  confirmatory hypotheses, a claim about one specific gene).
* `STAT-MT-004` **Use Benjamini-Yekutieli** when tests may be negatively or arbitrarily dependent
  [Benjamini 2001].
* `STAT-MT-005` **Independent filtering and weighting MUST use a covariate that is independent of
  the p-value under the null** (e.g. mean expression across all samples, not the fold change)
  [Ignatiadis 2016; Korthauer 2019].
* `STAT-MT-006` **Inspect the histogram of raw p-values** for every family. Expected shape: flat (true nulls) plus a
  peak near 0 (true effects). A hump in the middle or a U shape suggests a violated assumption or a wrong model. A
  peak near 1 also arises from valid but conservative tests: discrete tests, genes with very low counts, small
  samples. Before concluding that the model is miscalibrated, stratify the histogram by abundance (or another
  covariate) and, if needed, simulate under the fitted null model.
  * *Why:* when the model is correct and the test continuous, p-values of true nulls are uniform; deviations point
    to problems before any threshold is applied.
* `STAT-MT-007` **MUST NOT combine an FDR cut-off with an effect-size cut-off** (the "volcano plot" selection:
  adjusted $p < 0.05$ and $|\log_2 FC| > 1$) and claim FDR control for the result; the double
  filter inflates the FDR of the selected set [Ebrahimpoor 2021]. To require a minimum effect,
  test against the threshold (e.g. TREAT, $H_0: |\log_2 FC| \le \tau$) [McCarthy 2009].
* `STAT-MT-008` **Contrasts within one model** (pairwise comparisons, comparisons with a control) use
  procedures designed for them: Tukey for all pairwise, Dunnett for comparisons with a control,
  or a multivariate-t adjustment for arbitrary contrasts. These are less conservative than
  Bonferroni because they use the correlation between contrasts.
  * *Why:* Bonferroni is valid under any dependence (it uses the union bound) but does not exploit the known
    correlation between contrasts, so it is more conservative than necessary.
* `STAT-MT-009` Report adjusted values as "BH-adjusted p" or "q-value", stating the method; they are not
  interchangeable names for the same number unless the method is stated.

## 8. Effect sizes and uncertainty

* `STAT-EFF-001` **Every result reports an estimate with its uncertainty**: the effect size (difference in
  means, log fold change, odds ratio, correlation) and a confidence interval, plus the exact
  p-value when a test was done [Nakagawa 2007; Ho 2019].
  * *Why:* a p-value mixes the size of the effect with the sample size; the same p-value can correspond to a large
    effect in a small study or a trivial effect in a large one.
* `STAT-EFF-002` **p-values are graded evidence, not a verdict.** A p-value measures how incompatible the data
  are with the null model; it is not the probability that the null is true, and it does not
  measure the size or importance of an effect [Wasserstein 2016; Greenland 2016]. Thresholds
  (e.g. FDR < 0.05 to build a gene list) MAY be used when a decision is needed, but the
  underlying estimates are reported too.
* `STAT-EFF-003` **MUST NOT conclude "no effect" from $p > 0.05$.** Report the confidence interval; if it
  includes effects that would matter, the result is inconclusive. To claim the absence of a
  relevant effect, use an equivalence test against a pre-specified margin [Altman 1995;
  Lakens 2017].
* `STAT-EFF-004` **With large $n$, tiny p-values are expected for trivial effects** (thousands of cells, large
  cohorts). Judge results by effect size and its relevance, not by the p-value
  (Practitioner sources, Tang).
* `STAT-EFF-005` **Interpret the whole interval.** Values near its ends are less compatible with the data than
  values near the centre, but all are reasonably compatible [Amrhein 2019; Rafi 2020].
* `STAT-EFF-006` **Shrink noisy effect estimates for ranking and display** (e.g. log fold changes of low-count
  genes) with an established method, and say so [Zhu 2018; Stephens 2016].
  * *Trade-off:* shrunken estimates are biased towards zero; they are for ranking and display, while tests use the
    model's own estimates and standard errors.
* `STAT-EFF-007` **State which error bar is shown**: SD (standard deviation) describes the spread of the data; SEM (standard error of the mean) and confidence
  intervals describe the precision of the mean. Always give $n$ next to it [Cumming 2007].
  * *Why:* the SEM shrinks with $n$ while the SD does not, so the same data look more or less variable depending on
    the bar; an unlabelled bar cannot be interpreted.

## 9. Gene-set and functional enrichment

* `STAT-ENR-001` **The background (universe) MUST be the set of features that could have been selected**
  (e.g. genes expressed and tested in the experiment), not the whole genome [Wijesooriya 2022].
  *Why:* genes never measured can never be "significant"; including them makes any measured set
  look enriched.
* `STAT-ENR-002` **Correct for multiple testing across gene sets** (§7, multiple testing) and report the number
  of sets tested.
* `STAT-ENR-003` **Be aware of detection biases**: long or highly expressed genes are more often detected as
  DE; methods that correct for them exist and SHOULD be used when the bias is plausible
  [Timmons 2015].
* `STAT-ENR-004` **State the null hypothesis of the test** [Goeman 2007]:
  * *Competitive*: the genes in the set are more associated with the condition than other genes.
  * *Self-contained*: at least one gene in the set is associated with the condition.
* `STAT-ENR-005` Report the tool, its version, the gene-set database and its release, and the background list.

## 10. Resampling and simulation

* `STAT-RES-001` **Permutation tests MUST permute what is exchangeable under the null**: the condition labels
  of experimental units, within strata (batch, pair) when the design has them. Permuting cells
  instead of samples reproduces pseudoreplication (§4, hierarchical data).
  * *Why:* permuting units that are not exchangeable under the null destroys structure that the null still
    contains, producing a null distribution that is too narrow and p-values that are too small.
* `STAT-RES-002` **Exchangeability is a property of the design, not of the procedure.** Permuting within strata keeps the strata;
  it does not keep the association of the labels with other covariates (age, sex, RNA quality). When such
  covariates are associated with the labels, adding them to the model does not make plain label permutation valid;
  use restricted permutations (e.g. within matched sets), a method designed for nuisance covariates in that model
  class, or simulation from the fitted null model, and state the null that was constructed.
* `STAT-RES-003` **Permutation p-values are $p = (b + 1)/(m + 1)$**, where $b$ is the number of permutations with
  a statistic at least as extreme as the observed one and $m$ the number of permutations; they
  are never 0 [Phipson 2010]. "As extreme" needs a definition: for two-sided tests use twice the smaller one-sided
  p-value (capped at 1), as SciPy does; comparing absolute values is valid only for statistics centred at 0 under
  the null.
  * Recipe: `recipes/stats/permutation.py` (permutation within strata and this p-value).
* `STAT-RES-004` **The number of permutations limits the smallest attainable p-value** ($1/(m+1)$). Choose $m$
  so that this minimum is well below the threshold after multiple-testing correction. When the permutation result
  is a frequency (e.g. how often a pipeline makes any discovery), report its Monte Carlo uncertainty: with 0 events
  in 20 permutations the one-sided 95% upper bound is about 14%.
* `STAT-RES-005` **Bootstrap MUST resample experimental units**, not measurements within them. Prefer BCa (bias-corrected and accelerated) or
  percentile intervals from an established implementation.
  * *Why:* resampling measurements within a unit ignores the variation between units, so the intervals are too
    narrow.
* `STAT-RES-006` **Simulation studies state their data-generating model** and include the null scenario, to
  check that false-positive rates match the nominal level.
  * *Why:* a method that reports effects where none exist is invalid whatever its power; the null scenario is the
    first check.
* `STAT-RES-007` Seeds and determinism follow `CODING_GUIDELINES.md` §15 (reproducibility).

## 11. Bayesian analysis (optional)

* `STAT-BAYES-001` Bayesian models MAY be used when they fit the problem better (hierarchical shrinkage, prior
  information, small samples with strong structure).
* `STAT-BAYES-002` When used, state the priors and their justification, check convergence (R-hat, effective sample
  size), run posterior predictive checks, and report posterior intervals with their probability
  level [van de Schoot 2021].
* `STAT-BAYES-003` MUST NOT mix frequentist and Bayesian vocabulary: a credible interval is not a confidence
  interval, and a posterior probability is not a p-value.

## 12. What to report

Every statistical result reported in a notebook, table, figure or text MUST make the following
recoverable, in the result itself or in a linked methods note [Lang 2015; Percie du Sert 2020]:

1. `STAT-REP-001` The question as a comparison and the experimental unit.
2. `STAT-REP-002` $n$ per group (experimental units), and units excluded with the reason.
3. `STAT-REP-003` The model or test, with its formula (`~ batch + condition`), the software and its version.
4. `STAT-REP-004` One- or two-sided, and the multiple-testing method with the family size $m$.
5. `STAT-REP-005` The estimate, its confidence interval and the exact p-value (or adjusted p); not only stars
   or "p < 0.05".
6. `STAT-REP-006` The universe of the result, as required by `CODING_GUIDELINES.md` §11 (numerical code): input
   set, filters and parameters.
7. `STAT-REP-007` For figures: what each point, bar and error bar represents (`FIGURE_GUIDELINES.md`).

## 13. Implementation in Python

* `STAT-IMPL-001` Use established implementations; MUST NOT reimplement tests, models or corrections by hand
  except for teaching or validation.

| Task | Python | Notes |
|---|---|---|
| Linear / generalized linear models | `statsmodels.formula.api` (`ols`, `glm`) | Robust errors: `.fit(cov_type="HC3")`. |
| Negative binomial with estimated dispersion | `statsmodels.formula.api.negativebinomial` | `glm` with `NegativeBinomial()` family uses a fixed dispersion. |
| Mixed models (linear) | `statsmodels.formula.api.mixedlm` | For GLMMs, use R (`lme4`, `glmmTMB`) through a file interface. |
| GEE | `statsmodels.formula.api.gee` | State the working correlation structure. |
| Ordinal outcomes | `statsmodels.miscmodels.ordinal_model.OrderedModel` | |
| Welch t-test | `scipy.stats.ttest_ind(a, b, equal_var=False)` | |
| Rank tests | `scipy.stats.mannwhitneyu`, `scipy.stats.kruskal` | Interpret as in §3.3 (rank tests). |
| Tukey, Dunnett | `scipy.stats.tukey_hsd`, `scipy.stats.dunnett` | |
| Equivalence (TOST, two one-sided tests) | `statsmodels.stats.weightstats.ttost_ind` | |
| Multiple testing (BH, BY: Benjamini-Yekutieli, Holm) | `statsmodels.stats.multitest.multipletests` | `method="fdr_bh"`, `"fdr_by"`, `"holm"`. |
| Permutation tests, bootstrap | `scipy.stats.permutation_test`, `scipy.stats.bootstrap` | Randomized permutations already use the formula of §10 (resampling); custom loops MUST apply it. |
| Power and sample size | `statsmodels.stats.power` | |
| Bulk RNA-seq DE | `pydeseq2` | See `RNASEQ_BULK_GUIDELINES.md`; validate against DESeq2 (`CODING_GUIDELINES.md` §3, languages). |

* `STAT-IMPL-002` When no adequate Python implementation exists (IHW: independent hypothesis weighting; `qvalue`; `limma`; `emmeans`; GLMMs: generalized linear mixed models), call
  the R reference implementation through a standalone script, as `CODING_GUIDELINES.md` §3
  (languages) prescribes.
* `STAT-IMPL-003` Analysis parameters (thresholds, $\alpha$, number of permutations, seeds) are named constants in
  the configuration module (`CODING_GUIDELINES.md` §11, numerical code).

## 14. Glossary

* **Experimental unit**: the entity independently assigned to a condition; defines $n$.
* **Null hypothesis ($H_0$)**: the model of "no effect" against which the data are compared.
* **p-value**: probability, if $H_0$ and all model assumptions were true, of a result at least as
  extreme as the one observed.
* **Confidence interval (95%)**: range produced by a procedure that contains the true value in
  95% of repeated experiments; the values reasonably compatible with the data under the model.
* **Type I error / false positive**: declaring an effect that does not exist. **Type II error /
  false negative**: missing a real effect. **Power**: $1 -$ probability of a type II error.
* **FWER (family-wise error rate)**: probability of at least one false positive in the family. **FDR**: expected
  proportion of false positives among the results declared significant.
* **Effect size**: the magnitude of the difference or association, in interpretable units.
* **Overdispersion**: variance larger than the Poisson model predicts; typical of biological
  counts.
* **Pseudoreplication**: analysing non-independent measurements as independent.
* **Double dipping (circular analysis)**: using the same data to select and to test.
* **Data leakage**: information from the evaluation data reaching model training.
