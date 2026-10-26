# Statistics guidelines: references

> Companion to `STATISTICS_GUIDELINES.md` in <https://github.com/alexmascension/agent-coding-guidelines>.
> Compiled 2026-10-13. DOIs and bibliographic data verified against Crossref; BibTeX in `STATISTICS_REFERENCES.bib`.

Bracketed keys in the guides, e.g. `[Squair 2021]` in `STATISTICS_GUIDELINES.md` or `[Squair 2021†]` in the omics
guides, refer to the bold key at the start of each entry. Each entry gives the reference and the practical rule it supports. Entries are grouped by topic; within a topic, the most
general sources come first.

## Peer-reviewed literature

### A. P-values, significance and interpretation

- **[Wasserstein 2016]** Wasserstein & Lazar (2016). The ASA Statement on p-Values: Context, Process, and Purpose. *The American Statistician*. [doi:10.1080/00031305.2016.1154108](https://doi.org/10.1080/00031305.2016.1154108)  
  A p-value measures incompatibility of data with a specified model; it is not the probability that the hypothesis is true, and it does not measure effect size or importance.
- **[Wasserstein 2019]** Wasserstein et al. (2019). Moving to a World Beyond “p < 0.05”. *The American Statistician*. [doi:10.1080/00031305.2019.1583913](https://doi.org/10.1080/00031305.2019.1583913)  
  Do not use 'statistically significant' as a binary verdict; report estimates, uncertainty and context.
- **[McShane 2019]** McShane et al. (2019). Abandon Statistical Significance. *The American Statistician*. [doi:10.1080/00031305.2018.1527253](https://doi.org/10.1080/00031305.2018.1527253)  
  Argues for treating p-values as one continuous piece of evidence among others (prior evidence, design, plausibility), not as a threshold rule.
- **[Amrhein 2019]** Amrhein et al. (2019). Scientists rise up against statistical significance. *Nature*. [doi:10.1038/d41586-019-00857-9](https://doi.org/10.1038/d41586-019-00857-9)  
  Non-significance does not show absence of effect; interpret the whole confidence (compatibility) interval.
- **[Greenland 2016]** Greenland et al. (2016). Statistical tests, P values, confidence intervals, and power: a guide to misinterpretations. *European Journal of Epidemiology*. [doi:10.1007/s10654-016-0149-3](https://doi.org/10.1007/s10654-016-0149-3)  
  Catalogue of 25 common misinterpretations of p-values, confidence intervals and power; reference for wording of results.
- **[Rafi 2020]** Rafi & Greenland (2020). Semantic and cognitive tools to aid statistical science: replace confidence and significance by compatibility and surprise. *BMC Medical Research Methodology*. [doi:10.1186/s12874-020-01105-9](https://doi.org/10.1186/s12874-020-01105-9)  
  Suggests 'compatibility interval' wording and S-values ($s = -\log_2 p$) to express evidence against a model.
- **[Benjamin 2017]** Benjamin et al. (2017). Redefine statistical significance. *Nature Human Behaviour*. [doi:10.1038/s41562-017-0189-z](https://doi.org/10.1038/s41562-017-0189-z)  
  Proposes $\alpha = 0.005$ for new discoveries; useful as context when choosing thresholds, not as a universal rule.
- **[Altman 2017]** Altman & Krzywinski (2017). Interpreting P values. *Nature Methods*. [doi:10.1038/nmeth.4210](https://doi.org/10.1038/nmeth.4210)  
  Short Points of Significance primer on what a p-value does and does not mean.
- **[Halsey 2015]** Halsey et al. (2015). The fickle P value generates irreproducible results. *Nature Methods*. [doi:10.1038/nmeth.3288](https://doi.org/10.1038/nmeth.3288)  
  p-values vary widely between replications of the same experiment; a single p-value is a weak basis for decisions.
- **[Gelman & Stern 2006]** Gelman & Stern (2006). The Difference Between “Significant” and “Not Significant” is not Itself Statistically Significant. *The American Statistician*. [doi:10.1198/000313006X152649](https://doi.org/10.1198/000313006X152649)  
  Comparing a significant with a non-significant result is not a test of their difference; test the difference directly.
- **[Nieuwenhuis 2011]** Nieuwenhuis et al. (2011). Erroneous analyses of interactions in neuroscience: a problem of significance. *Nature Neuroscience*. [doi:10.1038/nn.2886](https://doi.org/10.1038/nn.2886)  
  Half of reviewed neuroscience papers claimed an interaction without testing it; test interactions explicitly.
- **[Altman 1995]** Altman & Bland (1995). Statistics notes: Absence of evidence is not evidence of absence. *BMJ*. [doi:10.1136/bmj.311.7003.485](https://doi.org/10.1136/bmj.311.7003.485)  
  Absence of evidence is not evidence of absence; use confidence intervals or equivalence tests to claim 'no effect'.
- **[Lakens 2017]** Lakens (2017). Equivalence Tests. *Social Psychological and Personality Science*. [doi:10.1177/1948550617697177](https://doi.org/10.1177/1948550617697177)  
  Equivalence testing (TOST) to support a claim that an effect is smaller than a pre-specified bound.

### B. Design, replication, power and sample size

- **[Lazic 2018]** Lazic et al. (2018). What exactly is ‘N’ in cell culture and animal experiments?. *PLOS Biology*. [doi:10.1371/journal.pbio.2005282](https://doi.org/10.1371/journal.pbio.2005282)  
  Define the experimental unit (the entity independently assigned to a condition); N is the number of experimental units, not of measurements.
- **[Lazic 2010]** Lazic (2010). The problem of pseudoreplication in neuroscientific studies: is it affecting your analysis?. *BMC Neuroscience*. [doi:10.1186/1471-2202-11-5](https://doi.org/10.1186/1471-2202-11-5)  
  Pseudoreplication (treating dependent measurements as independent) inflates false positives; common in lab biology.
- **[Aarts 2014]** Aarts et al. (2014). A solution to dependency: using multilevel analysis to accommodate nested data. *Nature Neuroscience*. [doi:10.1038/nn.3648](https://doi.org/10.1038/nn.3648)  
  Nested data (cells within animals) need multilevel models; ignoring the hierarchy inflates type I error.
- **[Blainey 2014]** Blainey et al. (2014). Replication. *Nature Methods*. [doi:10.1038/nmeth.3091](https://doi.org/10.1038/nmeth.3091)  
  Distinguishes technical from biological replicates and what each type supports in inference.
- **[Button 2013]** Button et al. (2013). Power failure: why small sample size undermines the reliability of neuroscience. *Nature Reviews Neuroscience*. [doi:10.1038/nrn3475](https://doi.org/10.1038/nrn3475)  
  Low power lowers the probability that a significant result is true and inflates effect-size estimates.
- **[Ioannidis 2005]** Ioannidis (2005). Why Most Published Research Findings Are False. *PLoS Medicine*. [doi:10.1371/journal.pmed.0020124](https://doi.org/10.1371/journal.pmed.0020124)  
  Small studies, small effects, flexible analyses and many tested relationships reduce the positive predictive value of findings.
- **[Simmons 2011]** Simmons et al. (2011). False-Positive Psychology. *Psychological Science*. [doi:10.1177/0956797611417632](https://doi.org/10.1177/0956797611417632)  
  Undisclosed analytic flexibility (optional stopping, covariate choice, outlier rules) inflates false positives; pre-specify and disclose.
- **[Lakens 2022]** Lakens (2022). Sample Size Justification. *Collabra: Psychology*. [doi:10.1525/collabra.33267](https://doi.org/10.1525/collabra.33267)  
  Six ways to justify sample size (resource limits, a-priori power, precision, heuristics…); state which one was used.
- **[Gelman & Carlin 2014]** Gelman & Carlin (2014). Beyond Power Calculations. *Perspectives on Psychological Science*. [doi:10.1177/1745691614551642](https://doi.org/10.1177/1745691614551642)  
  Type S (sign) and type M (magnitude) errors: in low-power studies significant estimates are exaggerated and may have the wrong sign.
- **[Hoenig 2001]** Hoenig & Heisey (2001). The Abuse of Power. *The American Statistician*. [doi:10.1198/000313001300339897](https://doi.org/10.1198/000313001300339897)  
  Post-hoc ('observed') power is a transformation of the p-value and adds no information; do not report it.
- **[Schurch 2016]** Schurch et al. (2016). How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use?. *RNA*. [doi:10.1261/rna.053959.115](https://doi.org/10.1261/rna.053959.115)  
  48-replicate yeast RNA-seq: at least 6 biological replicates per condition, 12 to detect most DE genes; tool choice depends on replicate number.
- **[Gierliński 2015]** Gierliński et al. (2015). Statistical models for RNA-seq data derived from a two-condition 48-replicate experiment. *Bioinformatics*. [doi:10.1093/bioinformatics/btv425](https://doi.org/10.1093/bioinformatics/btv425)  
  Same 48-replicate experiment: gene counts are broadly consistent with negative binomial (and log-normal) models, but a few bad replicates distort results and must be identified.
- **[Percie du Sert 2020]** Percie du Sert et al. (2020). The ARRIVE guidelines 2.0: Updated guidelines for reporting animal research. *PLOS Biology*. [doi:10.1371/journal.pbio.3000410](https://doi.org/10.1371/journal.pbio.3000410)  
  ARRIVE 2.0 'essential 10': design, sample size, inclusion/exclusion, randomisation, blinding, outcome measures, statistical methods.

### C. Choice of tests and models

- **[Delacre 2017]** Delacre et al. (2017). Why Psychologists Should by Default Use Welch’s t-test Instead of Student’s t-test. *International Review of Social Psychology*. [doi:10.5334/irsp.82](https://doi.org/10.5334/irsp.82)  
  Use Welch's t-test by default; Student's t-test is not robust to unequal variances with unequal group sizes.
- **[Rochon 2012]** Rochon et al. (2012). To test or not to test: Preliminary assessment of normality when comparing two independent samples. *BMC Medical Research Methodology*. [doi:10.1186/1471-2288-12-81](https://doi.org/10.1186/1471-2288-12-81)  
  Pre-testing normality before choosing a t-test or a rank test distorts error rates; choose the method a priori.
- **[Fagerland 2012]** Fagerland (2012). t-tests, non-parametric tests, and large studies—a paradox of statistical practice?. *BMC Medical Research Methodology*. [doi:10.1186/1471-2288-12-78](https://doi.org/10.1186/1471-2288-12-78)  
  With large samples the t-test is valid for non-normal data (CLT), while rank tests answer a different question.
- **[Lumley 2002]** Lumley et al. (2002). The Importance of the Normality Assumption in Large Public Health Data Sets. *Annual Review of Public Health*. [doi:10.1146/annurev.publhealth.23.100901.140546](https://doi.org/10.1146/annurev.publhealth.23.100901.140546)  
  In large samples linear regression and t-tests do not require normally distributed outcomes.
- **[Divine 2018]** Divine et al. (2018). The Wilcoxon–Mann–Whitney Procedure Fails as a Test of Medians. *The American Statistician*. [doi:10.1080/00031305.2017.1305291](https://doi.org/10.1080/00031305.2017.1305291)  
  Wilcoxon-Mann-Whitney does not test medians except under restrictive assumptions (shift model); it tests stochastic dominance.
- **[Krzywinski 2014a]** Krzywinski & Altman (2014). Nonparametric tests. *Nature Methods*. [doi:10.1038/nmeth.2937](https://doi.org/10.1038/nmeth.2937)  
  Points of Significance primer on rank-based tests, their assumptions and power.
- **[Altman 2006]** Altman & Royston (2006). The cost of dichotomising continuous variables. *BMJ*. [doi:10.1136/bmj.332.7549.1080](https://doi.org/10.1136/bmj.332.7549.1080)  
  Dichotomising continuous variables loses information and power and can create spurious effects.
- **[Vickers 2001]** Vickers & Altman (2001). Analysing controlled trials with baseline and follow up measurements. *BMJ*. [doi:10.1136/bmj.323.7321.1123](https://doi.org/10.1136/bmj.323.7321.1123)  
  With baseline and follow-up measurements, use ANCOVA (adjust for baseline) rather than change scores or comparing p-values.
- **[Bates 2015]** Bates et al. (2015). Fitting Linear Mixed-Effects Models Using lme4. *Journal of Statistical Software*. [doi:10.18637/jss.v067.i01](https://doi.org/10.18637/jss.v067.i01)  
  lme4: reference implementation of linear mixed models (random effects for nested/crossed designs).
- **[Barr 2013]** Barr et al. (2013). Random effects structure for confirmatory hypothesis testing: Keep it maximal. *Journal of Memory and Language*. [doi:10.1016/j.jml.2012.11.001](https://doi.org/10.1016/j.jml.2012.11.001)  
  For confirmatory tests, include the random-effects structure justified by the design; too-simple structures inflate type I error.
- **[Li 2022]** Li et al. (2022). Exaggerated false positives by popular differential expression methods when analyzing human population samples. *Genome Biology*. [doi:10.1186/s13059-022-02648-4](https://doi.org/10.1186/s13059-022-02648-4)  
  Permutation analysis: DESeq2/edgeR can exceed nominal FDR in large human population cohorts; authors recommend Wilcoxon for large n. Contested; treat as a warning to check FDR calibration (e.g. with permuted labels) in large heterogeneous cohorts.

### D. Pseudoreplication and single-cell differential expression

- **[Squair 2021]** Squair et al. (2021). Confronting false discoveries in single-cell differential expression. *Nature Communications*. [doi:10.1038/s41467-021-25960-2](https://doi.org/10.1038/s41467-021-25960-2)  
  Single-cell DE methods that ignore variation between biological replicates produce hundreds of false DE genes; pseudobulk methods perform best.
- **[Zimmerman 2021]** Zimmerman et al. (2021). A practical solution to pseudoreplication bias in single-cell studies. *Nature Communications*. [doi:10.1038/s41467-021-21038-1](https://doi.org/10.1038/s41467-021-21038-1)  
  Mixed models with a random effect for individual correct pseudoreplication in single-cell data.
- **[Crowell 2020]** Crowell et al. (2020). muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. *Nature Communications*. [doi:10.1038/s41467-020-19894-4](https://doi.org/10.1038/s41467-020-19894-4)  
  muscat: framework and benchmark for multi-sample, multi-condition single-cell DE; aggregation (pseudobulk) is competitive and fast.
- **[Junttila 2022]** Junttila et al. (2022). Benchmarking methods for detecting differential states between conditions from multi-subject single-cell RNA-seq data. *Briefings in Bioinformatics*. [doi:10.1093/bib/bbac286](https://doi.org/10.1093/bib/bbac286)  
  Benchmark of differential-state methods: naive cell-level methods suffer pseudoreplication bias; pseudobulk and mixed models are preferred.
- **[Murphy 2023]** Murphy et al. (2023). Avoiding false discoveries in single-cell RNA-seq by revisiting the first Alzheimer’s disease dataset. *eLife*. [doi:10.7554/eLife.90214](https://doi.org/10.7554/eLife.90214)  
  Re-analysis of a landmark Alzheimer snRNA-seq study: correct QC and pseudobulk DE change the reported DE genes substantially.
- **[Soneson 2018]** Soneson & Robinson (2018). Bias, robustness and scalability in single-cell differential expression analysis. *Nature Methods*. [doi:10.1038/nmeth.4612](https://doi.org/10.1038/nmeth.4612)  
  Benchmark of single-cell DE methods: results depend strongly on method and filtering; methods developed for bulk perform comparably or better.
- **[Ahlmann-Eltze 2023]** Ahlmann-Eltze & Huber (2023). Comparison of transformations for single-cell RNA-seq data. *Nature Methods*. [doi:10.1038/s41592-023-01814-1](https://doi.org/10.1038/s41592-023-01814-1)  
  Comparison of single-cell transformations: a shifted log with a sensible pseudocount followed by PCA performs as well as more complex methods.

### E. Batch effects

- **[Leek 2010]** Leek et al. (2010). Tackling the widespread and critical impact of batch effects in high-throughput data. *Nature Reviews Genetics*. [doi:10.1038/nrg2825](https://doi.org/10.1038/nrg2825)  
  Batch effects are widespread; design so that batch is not confounded with the variable of interest, and record processing metadata.
- **[Nygaard 2015]** Nygaard et al. (2015). Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses. *Biostatistics*. [doi:10.1093/biostatistics/kxv027](https://doi.org/10.1093/biostatistics/kxv027)  
  Batch correction 'preserving group differences' on unbalanced designs yields exaggerated confidence; include batch in the model instead.
- **[Zhang 2020]** Zhang et al. (2020). ComBat-seq: batch effect adjustment for RNA-seq count data. *NAR Genomics and Bioinformatics*. [doi:10.1093/nargab/lqaa078](https://doi.org/10.1093/nargab/lqaa078)  
  ComBat-seq: batch adjustment on counts; corrected data are for visualisation or for methods that cannot model batch.

### F. Circular analysis, double dipping and data leakage

- **[Kriegeskorte 2009]** Kriegeskorte et al. (2009). Circular analysis in systems neuroscience: the dangers of double dipping. *Nature Neuroscience*. [doi:10.1038/nn.2303](https://doi.org/10.1038/nn.2303)  
  Circular analysis: selecting data by a criterion and then testing that same criterion on the same data inflates results.
- **[Gao 2022]** Gao et al. (2022). Selective Inference for Hierarchical Clustering. *Journal of the American Statistical Association*. [doi:10.1080/01621459.2022.2116331](https://doi.org/10.1080/01621459.2022.2116331)  
  Testing for differences between clusters estimated on the same data is invalid with standard tests; selective inference corrects it.
- **[Neufeld 2022]** Neufeld et al. (2022). Inference after latent variable estimation for single-cell RNA sequencing data. *Biostatistics*. [doi:10.1093/biostatistics/kxac047](https://doi.org/10.1093/biostatistics/kxac047)  
  Count splitting: split counts into independent train/test parts to estimate latent structure and test associations validly.
- **[Lee 2025]** Lee et al. (2025). ClusterDE: A Statistical Software Package for Removing Double-Dipping Bias in Post-Clustering Differential Expression Analysis. *Journal of Computational Biology*. [doi:10.1177/15578666251383562](https://doi.org/10.1177/15578666251383562)  
  ClusterDE: synthetic null to remove double-dipping false positives in post-clustering DE.
- **[Ambroise 2002]** Ambroise & McLachlan (2002). Selection bias in gene extraction on the basis of microarray gene-expression data. *Proceedings of the National Academy of Sciences*. [doi:10.1073/pnas.102102699](https://doi.org/10.1073/pnas.102102699)  
  Feature selection outside cross-validation produces optimistic error rates; selection must be repeated inside each fold.
- **[Varma 2006]** Varma & Simon (2006). Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics*. [doi:10.1186/1471-2105-7-91](https://doi.org/10.1186/1471-2105-7-91)  
  Tuning hyperparameters with the same cross-validation used for error estimation is biased; use nested cross-validation.
- **[Kapoor 2023]** Kapoor & Narayanan (2023). Leakage and the reproducibility crisis in machine-learning-based science. *Patterns*. [doi:10.1016/j.patter.2023.100804](https://doi.org/10.1016/j.patter.2023.100804)  
  Taxonomy of data leakage in ML-based science (preprocessing on full data, duplicates, temporal leakage) and a model-info sheet to prevent it.
- **[Chari 2023]** Chari & Pachter (2023). The specious art of single-cell genomics. *PLOS Computational Biology*. [doi:10.1371/journal.pcbi.1011288](https://doi.org/10.1371/journal.pcbi.1011288)  
  t-SNE/UMAP embeddings distort distances and global structure; do not draw quantitative conclusions from 2-D embeddings.

### G. Multiple testing and false discovery rate

- **[Benjamini 1995]** Benjamini & Hochberg (1995). Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. *Journal of the Royal Statistical Society Series B: Statistical Methodology*. [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)  
  Benjamini-Hochberg procedure for FDR control; default for genome-scale testing.
- **[Benjamini 2001]** Benjamini & Yekutieli (2001). The control of the false discovery rate in multiple testing under dependency. *The Annals of Statistics*. [doi:10.1214/aos/1013699998](https://doi.org/10.1214/aos/1013699998)  
  Benjamini-Yekutieli: BH holds under positive dependence; BY controls FDR under arbitrary dependence at a cost in power.
- **[Storey 2003]** Storey & Tibshirani (2003). Statistical significance for genomewide studies. *Proceedings of the National Academy of Sciences*. [doi:10.1073/pnas.1530509100](https://doi.org/10.1073/pnas.1530509100)  
  q-values and estimation of the proportion of true nulls $\pi_0$ for genome-wide studies.
- **[Stephens 2016]** Stephens (2016). False discovery rates: a new deal. *Biostatistics*. [doi:10.1093/biostatistics/kxw041](https://doi.org/10.1093/biostatistics/kxw041)  
  Adaptive shrinkage (ashr): FDR and effect estimates using effect sizes and their standard errors, local false sign rate.
- **[Ignatiadis 2016]** Ignatiadis et al. (2016). Data-driven hypothesis weighting increases detection power in genome-scale multiple testing. *Nature Methods*. [doi:10.1038/nmeth.3885](https://doi.org/10.1038/nmeth.3885)  
  IHW: weight hypotheses using an independent covariate (e.g. mean expression) to gain power while controlling FDR.
- **[Korthauer 2019]** Korthauer et al. (2019). A practical guide to methods controlling false discoveries in computational biology. *Genome Biology*. [doi:10.1186/s13059-019-1716-1](https://doi.org/10.1186/s13059-019-1716-1)  
  Benchmark of modern FDR methods with covariates; informative covariates improve power; choose covariates independent of p-values under the null.
- **[Goeman 2014]** Goeman & Solari (2014). Multiple hypothesis testing in genomics. *Statistics in Medicine*. [doi:10.1002/sim.6082](https://doi.org/10.1002/sim.6082)  
  Review of multiple testing in genomics: FWER vs FDR, dependence, and what each error rate guarantees.
- **[Ebrahimpoor 2021]** Ebrahimpoor & Goeman (2021). Inflated false discovery rate due to volcano plots: problem and solutions. *Briefings in Bioinformatics*. [doi:10.1093/bib/bbab053](https://doi.org/10.1093/bib/bbab053)  
  Volcano-plot double filtering (adjusted p plus fold-change cutoff) inflates FDR in the selected subset; test against a threshold instead.
- **[McCarthy 2009]** McCarthy & Smyth (2009). Testing significance relative to a fold-change threshold is a TREAT. *Bioinformatics*. [doi:10.1093/bioinformatics/btp053](https://doi.org/10.1093/bioinformatics/btp053)  
  TREAT: test whether |log fold change| exceeds a threshold, instead of filtering significant genes by fold change post hoc.
- **[Phipson 2010]** Phipson & Smyth (2010). Permutation P-values Should Never Be Zero: Calculating Exact P-values When Permutations Are Randomly Drawn. *Statistical Applications in Genetics and Molecular Biology*. [doi:10.2202/1544-6115.1585](https://doi.org/10.2202/1544-6115.1585)  
  Permutation p-values must be $(b+1)/(m+1)$, never zero, where $b$ is the number of permutation statistics at least as extreme and $m$ the number of permutations.

### H. Effect sizes, uncertainty and their display

- **[Nakagawa 2007]** Nakagawa & Cuthill (2007). Effect size, confidence interval and statistical significance: a practical guide for biologists. *Biological Reviews*. [doi:10.1111/j.1469-185X.2007.00027.x](https://doi.org/10.1111/j.1469-185X.2007.00027.x)  
  Report effect sizes with confidence intervals; guide to effect-size measures for biologists.
- **[Ho 2019]** Ho et al. (2019). Moving beyond P values: data analysis with estimation graphics. *Nature Methods*. [doi:10.1038/s41592-019-0470-3](https://doi.org/10.1038/s41592-019-0470-3)  
  Estimation statistics: report the effect size and its bootstrap CI and show it in the figure (estimation plots).
- **[Krzywinski 2013a]** Krzywinski & Altman (2013). Importance of being uncertain. *Nature Methods*. [doi:10.1038/nmeth.2613](https://doi.org/10.1038/nmeth.2613)  
  Points of Significance: sampling distributions, standard error and why uncertainty must accompany estimates.
- **[Krzywinski 2013b]** Krzywinski & Altman (2013). Error bars. *Nature Methods*. [doi:10.1038/nmeth.2659](https://doi.org/10.1038/nmeth.2659)  
  Error bars: SD describes spread, SEM and CI describe precision; always state which is shown.
- **[Cumming 2007]** Cumming et al. (2007). Error bars in experimental biology. *The Journal of Cell Biology*. [doi:10.1083/jcb.200611141](https://doi.org/10.1083/jcb.200611141)  
  Rules for error bars: state n and the type of bar; descriptive vs inferential bars; overlap rules are often misread.
- **[Weissgerber 2015]** Weissgerber et al. (2015). Beyond Bar and Line Graphs: Time for a New Data Presentation Paradigm. *PLOS Biology*. [doi:10.1371/journal.pbio.1002128](https://doi.org/10.1371/journal.pbio.1002128)  
  Bar graphs of continuous data hide distributions; show individual points, especially for small n.
- **[Krzywinski 2013c]** Krzywinski & Altman (2013). Power and sample size. *Nature Methods*. [doi:10.1038/nmeth.2738](https://doi.org/10.1038/nmeth.2738)  
  Points of Significance on power and sample size, with relation between effect size, variance, n and $\alpha$.
- **[Krzywinski 2014b]** Krzywinski & Altman (2014). Comparing samples—part I. *Nature Methods*. [doi:10.1038/nmeth.2858](https://doi.org/10.1038/nmeth.2858)  
  Points of Significance: comparing two samples (t-tests, pairing).
- **[Krzywinski 2014c]** Krzywinski & Altman (2014). Comparing samples—part II. *Nature Methods*. [doi:10.1038/nmeth.2900](https://doi.org/10.1038/nmeth.2900)  
  Points of Significance: multiple comparisons and their correction.
- **[Zhu 2018]** Zhu et al. (2018). Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics*. [doi:10.1093/bioinformatics/bty895](https://doi.org/10.1093/bioinformatics/bty895)  
  apeglm: shrink log fold changes for ranking and visualisation; unshrunk estimates of low-count genes are noisy.

### I. Reference implementations for count data

- **[Love 2014]** Love et al. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)  
  DESeq2: negative binomial GLM with dispersion and fold-change shrinkage; reference for bulk RNA-seq DE.
- **[Law 2014]** Law et al. (2014). voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biology*. [doi:10.1186/gb-2014-15-2-r29](https://doi.org/10.1186/gb-2014-15-2-r29)  
  voom: precision weights allow linear models (limma) on log-CPM RNA-seq data.
- **[Smyth 2004]** Smyth (2004). Linear Models and Empirical Bayes Methods for Assessing Differential Expression in Microarray Experiments. *Statistical Applications in Genetics and Molecular Biology*. [doi:10.2202/1544-6115.1027](https://doi.org/10.2202/1544-6115.1027)  
  limma: empirical Bayes moderation of variances, important with few replicates.

### J. Gene-set enrichment

- **[Goeman 2007]** Goeman & Bühlmann (2007). Analyzing gene expression data in terms of gene sets: methodological issues. *Bioinformatics*. [doi:10.1093/bioinformatics/btm051](https://doi.org/10.1093/bioinformatics/btm051)  
  Gene-set tests differ in their null hypothesis (competitive vs self-contained) and in the sampling unit (genes vs subjects); choose deliberately.
- **[Timmons 2015]** Timmons et al. (2015). Multiple sources of bias confound functional enrichment analysis of global -omics data. *Genome Biology*. [doi:10.1186/s13059-015-0761-7](https://doi.org/10.1186/s13059-015-0761-7)  
  Over-representation analysis is biased by an inappropriate background and by detection biases such as gene length and expression level.
- **[Wijesooriya 2022]** Wijesooriya et al. (2022). Urgent need for consistent standards in functional enrichment analysis. *PLOS Computational Biology*. [doi:10.1371/journal.pcbi.1009935](https://doi.org/10.1371/journal.pcbi.1009935)  
  Screen of 186 papers: 95% of over-representation analyses lacked an appropriate background list; many lacked FDR correction or methodological detail.

### K. Compositional data

- **[Gloor 2017]** Gloor et al. (2017). Microbiome Datasets Are Compositional: And This Is Not Optional. *Frontiers in Microbiology*. [doi:10.3389/fmicb.2017.02224](https://doi.org/10.3389/fmicb.2017.02224)  
  Relative-abundance data are compositional; analyse with log-ratio methods rather than correlations of proportions.

### L. General practice and reporting

- **[Kass 2016]** Kass et al. (2016). Ten Simple Rules for Effective Statistical Practice. *PLOS Computational Biology*. [doi:10.1371/journal.pcbi.1004961](https://doi.org/10.1371/journal.pcbi.1004961)  
  Ten rules: statistical methods should enable data to answer scientific questions; plan early; signals come with noise; check assumptions; be reproducible.
- **[Si 2025]** Si, Hashimoto & Yang (2025). The Ideation-Execution Gap: Execution Outcomes of LLM-Generated versus Human Research Ideas. *arXiv preprint* 2506.20803 (not peer reviewed; bibliographic data checked against arXiv on 2026-10-26). [doi:10.48550/arXiv.2506.20803](https://doi.org/10.48550/arXiv.2506.20803)  
  43 researchers executed randomly assigned LLM-generated or expert research ideas; after execution the review scores of the LLM ideas fell significantly more than those of the expert ideas on every metric, closing the gap seen at the idea stage. How promising an idea looks before execution is a poor guide to its outcome; cited in `AGENTIC_GUIDELINES.md` (idea triage).
- **[Makin 2019]** Makin & Orban de Xivry (2019). Ten common statistical mistakes to watch out for when writing or reviewing a manuscript. *eLife*. [doi:10.7554/eLife.48175](https://doi.org/10.7554/eLife.48175)  
  Ten common mistakes: absent controls, comparing effects without testing interaction, inflating units, spurious correlations, circularity, small samples, flexibility, no multiple-comparison correction, over-interpreting non-significance, correlation vs causation.
- **[Lang 2015]** Lang & Altman (2015). Basic statistical reporting for articles published in Biomedical Journals: The “Statistical Analyses and Methods in the Published Literature” or the SAMPL Guidelines. *International Journal of Nursing Studies*. [doi:10.1016/j.ijnurstu.2014.09.006](https://doi.org/10.1016/j.ijnurstu.2014.09.006)  
  SAMPL guidelines: what to report for each statistical method (descriptives, test, sidedness, exact p, CIs, software).
- **[Lever 2017]** Lever et al. (2017). Principal component analysis. *Nature Methods*. [doi:10.1038/nmeth.4346](https://doi.org/10.1038/nmeth.4346)  
  Points of Significance on PCA: centring and scaling choices, interpretation of components and explained variance.
- **[van de Schoot 2021]** van de Schoot et al. (2021). Bayesian statistics and modelling. *Nature Reviews Methods Primers*. [doi:10.1038/s43586-020-00001-2](https://doi.org/10.1038/s43586-020-00001-2)  
  Primer on Bayesian modelling workflow: priors, computation, model checking, reporting.
- **[Sandve 2013]** Sandve et al. (2013). Ten Simple Rules for Reproducible Computational Research. *PLoS Computational Biology*. [doi:10.1371/journal.pcbi.1003285](https://doi.org/10.1371/journal.pcbi.1003285)  
  Reproducible computational research: record versions, seeds, intermediate results and the exact commands.

## Practitioner sources

LinkedIn posts require login and could not be retrieved; the positions below come from the authors' public repositories
and blog, which cover the same topics as their posts.

### Adrian Olszewski (clinical-trials biostatistician)

Sources: [model-based-testing-hypotheses](https://github.com/adrianolszewski/model-based-testing-hypotheses),
[Logistic-regression-is-regression](https://github.com/adrianolszewski/Logistic-regression-is-regression),
[Mann-Whitney-is-not-about-medians-in-general](https://github.com/adrianolszewski/Mann-Whitney-is-not-about-medians-in-general),
and his LinkedIn *#StatSins* series (titles visible in search results).

- **Test hypotheses through models.** Most classical tests are special cases of (generalised) linear models. Fitting the
  model (GLM, GLS, GEE, mixed model) and testing on it allows covariates, multi-way interactions and planned contrasts
  that standalone tests cannot handle.
- **Fix assumption violations with the right model, not by switching to a rank test.** Heteroscedasticity → GLS or
  robust (sandwich) errors; dependence → GEE or mixed models; skewed positive outcomes → gamma GLM; counts → Poisson or
  negative binomial; proportions → binomial or beta; ordinal outcomes → proportional-odds model.
- **Wald, likelihood-ratio and score tests can disagree** in small samples or near boundaries. Marked disagreement
  signals that large-sample approximations are unreliable; the likelihood-ratio test is generally preferred there.
- **Contrasts after a model** should use multiplicity adjustments that exploit the joint distribution of the estimates
  (multivariate t), which are less conservative than Bonferroni.
- **Mann-Whitney, Kruskal-Wallis and Friedman do not compare medians** in general; they test stochastic dominance.
  Mann-Whitney is equivalent to a proportional-odds (ordinal logistic) model with one binary predictor.
- **Logistic regression is a regression** of the probability of the outcome; classification is a separate thresholding
  step. Use it (and its ordinal, multinomial, conditional and mixed variants) for inference, not only prediction.
- *#StatSins* themes: do not transform variables or delete outliers to force normality; do not rely blindly on the CLT;
  p-values and confidence intervals are routinely misinterpreted (he recommends Greenland et al. 2016, section A).

### Ming "Tommy" Tang (computational biologist, *Chatomics* blog)

Sources: posts on [divingintogeneticsandgenomics.com](https://divingintogeneticsandgenomics.com/) (listed below) and the
[RNA-seq-analysis](https://github.com/crazyhottommy/RNA-seq-analysis) and
[scRNAseq-analysis-notes](https://github.com/crazyhottommy/scRNAseq-analysis-notes) repositories.

- **Large N makes p-values uninformative**; with thousands of cells or data points, report and judge effect sizes
  (log fold change, correlation coefficient). ([Bioinformatics is not (just) statistics](https://divingintogeneticsandgenomics.com/post/bioinformatics-is-not-just-statistics/);
  [How to test if two distributions are different](https://divingintogeneticsandgenomics.com/post/how-to-test-if-two-distributions-are-different/))
- **Double dipping**: clustering and then testing markers between those clusters gives deflated p-values; treat them as
  ranking scores, not evidence. ([Do you really understand log2FC in scRNA-seq](https://divingintogeneticsandgenomics.com/post/do-you-really-understand-log2fold-change-in-single-cell-rnaseq-data/);
  [scRNA-seq clustering significance test](https://divingintogeneticsandgenomics.com/post/scrnaseq-clustering-significant-test-an-unsolvable-problem/))
- **Batch in the design, not in a pre-corrected matrix.** Include batch as a covariate in limma/DESeq2; use corrected
  matrices only for visualisation; check from metadata whether "replicate" really corresponds to a shared technical
  batch. ([Batch effect: correct or not](https://divingintogeneticsandgenomics.com/post/batch-effect-correct-or-not-bulk-rnaseq/))
- **Pseudobulk for multi-sample single-cell comparisons**; cells are not independent replicates.
  ([How to create pseudobulk](https://divingintogeneticsandgenomics.com/post/how-to-create-pseudobulk-from-single-cell-rnaseq-data/))
- **Multiple testing**: state the null hypothesis; inspect the p-value histogram; BH or Storey FDR for gene lists,
  Bonferroni/FWER where any false positive is costly; q-values are per-gene.
  ([Understanding p-value, multiple comparisons, FDR and q-value](https://divingintogeneticsandgenomics.com/post/understanding-p-value-multiple-comparisons-fdr-and-q-value/))
- **Check implementation details** (pseudocounts in log fold change, CLR variants, defaults) against a manual calculation
  for a few genes. ([CLR normalization details](https://divingintogeneticsandgenomics.com/post/details-in-centered-log-ratio-clr-normalization-for-cite-seq-protein-count-data/))
- **Zeros must be explicit**: groups with zero counts (e.g. a cell type absent in a sample) must appear as zero rows;
  test abundance on counts, not percentages. ([dplyr count misses factor levels](https://divingintogeneticsandgenomics.com/post/dplyr-count-misses-factor-levels-a-case-in-comparing-scrnaseq-cell-type-abundance/))
- **Missing is not zero**: distinguish undetected from unassayed values before analysis.
  ([PD-L1 RNA vs protein](https://divingintogeneticsandgenomics.com/post/are-pdl1-rna-and-protein-levels-correlated-in-cancer-cell-lines/))
- **Use permutation** to calibrate statistics without reliable asymptotic nulls (distribution distances, number of PCs).
  ([Permutation test for PCA components](https://divingintogeneticsandgenomics.com/post/permute-test-for-pca-components/))
- **High-dimensional prediction**: regularise when $p \gg n$; keep the test set untouched until final evaluation.
  ([Marker selection with regularised logistic regression](https://divingintogeneticsandgenomics.com/post/marker-gene-selection-using-logistic-regression-and-regularization-for-scrnaseq/))
- **Understand the design and run sanity checks** against known positive and negative controls before interpreting.
  ([Hidden skills beyond programming](https://divingintogeneticsandgenomics.com/post/hidden-skills-beyond-programming-for-computational-biology/))

## Books and online resources (not verified in this compilation)

- Harrell, F. E. *Regression Modeling Strategies* (2nd ed., Springer, 2015) and *Biostatistics for Biomedical Research* (online).
- Holmes, S. & Huber, W. *Modern Statistics for Modern Biology* (Cambridge University Press, 2019; free online).
- Irizarry, R. A. & Love, M. I. *Data Analysis for the Life Sciences* (Leanpub; free online).
- Lazic, S. E. *Experimental Design for Laboratory Biologists* (Cambridge University Press, 2016).
- Lindeløv, J. K. *Common statistical tests are linear models* (online cheat sheet).
- Krzywinski, M. & Altman, N. *Points of Significance* column, *Nature Methods* (2013–), several entries cited above.
