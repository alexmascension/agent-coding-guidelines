# Bulk RNA-seq guidelines

> Source: <https://github.com/alexmascension/agent-coding-guidelines> · Version 2.1.0 (2026-10-14)
>
> Key words **MUST**, **SHOULD**, **MAY** follow RFC 2119. Bracketed keys such as `[Conesa 2016]` point to
> `OMICS_REFERENCES.md`; keys with a dagger, e.g. `[Schurch 2016†]`, point to `STATISTICS_REFERENCES.md`. Each rule is
> followed by ***Why*** (mechanism), ***Trade-off*** (cost) and ***Deviate when*** (justified exceptions) where they
> apply.

## 0. How to use this document

* Audience: coding agents and people with introductory statistics.
* Copy with `OMICS_REFERENCES.md`, `STATISTICS_GUIDELINES.md` and `STATISTICS_REFERENCES.md`. Related guides:
  `SEX_INFERENCE_GUIDELINES.md` (sex checks and sex-matched references), `SAMPLE_IDENTITY_GUIDELINES.md` (duplicates,
  swaps), `DECONVOLUTION_GUIDELINES.md` (cell composition).
* Scope: short-read bulk RNA-seq of human and mouse, from FASTQ to gene- and transcript-level differential analysis and
  enrichment. Workflow references: [Conesa 2016; Love 2016a; Chen 2016; Law 2014†].
* Language: Python first (`pydeseq2`, `gseapy`, `decoupler`); R packages through scripts (`CODING_GUIDELINES.md` §3).

## 1. Design

* **The experimental unit is the biological sample**; technical replicates (re-sequencing of a library) are summed.
* **At least 6 biological replicates per group** for a general DE (differential expression) study [Schurch 2016†].
  * *Why:* with fewer, only large changes are detected, and estimates of dispersion are unstable.
* **More replicates rather than more depth** beyond about 20 million reads per sample [Ching 2014].
* **Balance technical factors across conditions** (extraction, library and sequencing batches) and record them with
  RIN (RNA integrity number) or DV200 (percentage of RNA fragments longer than 200 nucleotides), library kit, strandedness, read length, and, for postmortem tissue, postmortem interval.
  * *Why:* a technical factor confounded with the condition cannot be removed afterwards [Leek 2010†]; RNA quality is
    the main technical factor in postmortem tissue (§7).
* **Record sex, age and other biological covariates** and plan them in the model.

## 2. Reference genome and annotation

* **One genome build and one annotation release for all samples**, recorded in `PROJECT_CONTEXT.md` §6.
  * *Why:* the gene model changes counts substantially: comparing RefSeq and Ensembl annotations, only 16% of genes had
    identical counts and 9% differed by 50% or more [Zhao 2015]. Changing the annotation mid-project requires
    re-quantifying everything.
* **Sex-matched references** (Y masked for XX; Y PAR: pseudoautosomal region, masked for XY) are recommended when sex-chromosome or PAR genes
  matter; routing follows `SEX_INFERENCE_GUIDELINES.md` §4–6 [Olney 2020].
* **Full annotation, filtered at analysis time**: quantify against the complete transcript annotation and apply
  filters (e.g. primary or MANE transcripts, expression) afterwards.
  * *Why:* removing transcripts from the index forces their reads onto similar isoforms, biasing those isoforms.
  * *Trade-off:* more ambiguity between isoforms, larger outputs; inferential replicates (§4) quantify that ambiguity.
  * *Deviate when:* for transcript-usage analyses, prefiltering lowly used isoforms improves FDR (false discovery rate) control [Soneson 2016b];
    it is done on the expression estimates, after quantification.
* **Identify genes by stable Ensembl IDs**; map to symbols only for display, through the same release.

## 3. Read-level quality control

* FastQC per FASTQ, aggregated with MultiQC [Ewels 2016]; keep the report with the run outputs.
* Per sample: total reads, mapping rate, rRNA fraction, duplication, gene-body coverage (3' bias indicates RNA
  degradation) and **inferred strandedness**, which MUST match the library kit.
* Sex and identity checks follow `SEX_INFERENCE_GUIDELINES.md` and `SAMPLE_IDENTITY_GUIDELINES.md`; run them before
  any analysis that combines samples.

## 4. Quantification

* **Default: `nf-core/rnaseq` with `--aligner star_salmon`** [Ewels 2020], pinned to a release.
  * *What it does:* STAR aligns to the genome [Dobin 2012] and writes a second BAM projected onto the transcriptome;
    Salmon in alignment mode distributes each fragment among compatible transcripts with an EM (expectation-maximization) algorithm [Patro 2017].
  * *Why:* genome BAMs are needed for QC (quality control), sex and identity fingerprints, splice junctions and regions of interest,
    and the transcript-level estimates come from the same alignments.
  * *Deviate when:* no BAMs are needed: Salmon in selective-alignment mode with the genome as decoys is faster and
    nearly as accurate [Srivastava 2020]; kallisto is similar [Bray 2016]. In benchmarks the differences between
    modern quantifiers are smaller than those caused by the annotation [Teng 2016; Zhang 2017; Zhao 2015].
* **Salmon options: `--seqBias --gcBias` always; `--numGibbsSamples 20` (or bootstraps) whenever transcript-level
  analysis is possible.** With nf-core: `--extra_salmon_quant_args '--seqBias --gcBias --numGibbsSamples 20'`.
  * *Why `--seqBias`:* random-hexamer priming over-represents some fragment starts; Salmon learns the pattern per
    sample and corrects effective transcript lengths.
  * *Why `--gcBias`:* PCR amplifies GC-extreme fragments less efficiently, and the bias differs between batches; if
    two isoforms differ in a GC-rich exon, uncorrected bias changes which isoform appears dominant and produces false
    isoform switches between batches [Love 2016b].
  * *Why inferential replicates:* when isoforms share most exons, the EM gives one point estimate whose uncertainty is
    technical. Gibbs samples quantify it, and methods that use them control the FDR for ambiguous transcripts: swish
    [Zhu 2019] and edgeR with divided counts for transcript expression [Baldoni 2023] and usage [Baldoni 2025].
  * *Trade-off:* minutes per sample and tens of MB of output per sample; far cheaper than re-quantifying a cohort.
  * *Deviate when:* only gene-level DE will ever be done; inferential replicates are then optional.
* **Transcript to gene**: aggregate with tximport [Soneson 2016a].
  * For DESeq2 or edgeR in R: counts with average-transcript-length offsets.
  * **For `pydeseq2`: `countsFromAbundance = "lengthScaledTPM"`**, because `pydeseq2` (0.5) accepts sample size
    factors but not gene-by-sample offsets.
  * *Why:* gene counts that ignore changes in isoform length between samples bias gene-level DE; both options correct
    for it [Soneson 2016a].
  * Recipe: `recipes/bulk/counts_io.py` (nf-core gene tables) and `recipes/bulk/tximport_lengthscaled.R` (Salmon
    runs outside nf-core).
* **DE tools take counts, never TPM (transcripts per million), FPKM (fragments per kilobase per million mapped fragments) or normalized values**; TPM and RPKM (reads per kilobase per million mapped reads) are within-sample measures and MUST NOT
  be used to compare samples or protocols [Zhao 2020].

## 5. Filtering and normalization

* **Filter low-count genes** with a design-based rule, e.g. a minimum count (around 10) in at least as many samples as
  the smallest group.
  * *Why:* low-count genes cannot reach significance and only increase the multiple-testing burden. The rule uses group
    size, not the condition of individual samples, so it does not bias the test.
  * Recipe: `recipes/bulk/deseq2.py` (`filter_genes`).
* **Use the normalization of the DE method**: median-of-ratios (DESeq2, `pydeseq2`) [Anders 2010] or TMM (trimmed mean of M-values; edgeR)
  [Robinson 2010]. Both assume that most genes are not DE.
  * *Deviate when:* a global shift in expression is expected (transcriptional amplification, very different
    composition); use spike-ins or control genes (RUV, removal of unwanted variation) [Risso 2014; Evans 2017].
* **Gene length bias can be sample-specific** and is not removed by TMM, median-of-ratios or RPKM; it creates false
  enrichment of short-gene (ribosomal) and long-gene (extracellular matrix) sets [Mandelboum 2019]. Check it by
  plotting log fold changes against gene length; if present, use conditional quantile normalization (cqn) or EDASeq
  with length as a sample-specific covariate.
  * Recipe: `recipes/bulk/length_bias.py` (per-sample length slopes, their association with the condition, and
    log fold change by length decile).

## 6. Exploration

* **Transformation for PCA (principal component analysis), heatmaps and clustering: VST (variance-stabilizing transformation)** (`pydeseq2`: `dds.vst(use_design=False)`).
  * *Why VST:* it stabilizes the variance across the expression range (raw or simply log-transformed counts are
    dominated by noise at low counts); for high counts it approaches $\log_2$ of normalized counts, so it remains
    interpretable; it is available in Python and scales to thousands of samples.
  * *Trade-off:* `use_design=False` (equivalent to `blind=TRUE` in DESeq2) is right for QC; for downstream use in
    designs with large expected differences, `use_design=True` avoids shrinking real differences.
  * *Deviate when:* rlog (R only, slow beyond a few dozen samples) for small studies; log-CPM with a prior count
    (edgeR/limma) when working in the limma framework.
* **PCA** on the top ~500 variable genes, coloured by condition and by **every** technical variable; sample correlation
  heatmap.
* Samples that separate by a technical variable are a modelling problem (§7), not outliers to drop. Exclude samples only
  by criteria fixed in advance (`STATISTICS_GUIDELINES.md` §2).

## 7. Known and hidden factors

* **Known covariates first**: batch, platform, sex, age, and RNA quality (RIN), plus postmortem interval for postmortem
  tissue.
  * *Why:* including RIN as a covariate corrects most degradation-driven bias, provided RIN is not confounded with the
    condition [Gallego Romero 2014].
* **Hidden factors MAY be estimated and added as covariates**, with every biological variable of interest already in
  the model. Choice of method:

  | Method | Use when | Main risk |
  |---|---|---|
  | **qSVA** (quality surrogate variable analysis) [Jaffe 2017] | Postmortem brain, RNA quality correlated with diagnosis | Developed in cortex; transfer to other tissues (e.g. spinal cord) is untested |
  | **RUV** with control genes, replicate samples (RUVs), or RUV-III with pseudo-replicates (PRPS) for large multi-batch cohorts [Gagnon-Bartsch 2011; Risso 2014; Molania 2022] | Negative controls or replicate samples exist | Poor controls absorb biology |
  | **SVA** (surrogate variable analysis) / svaseq [Leek 2007; Leek 2014] | No controls, moderate sample size | Removes real biological signal not in the model (e.g. sex, copy number) [Jaffe 2015] |
  | **PCA of expression** | QTL (quantitative trait locus) mapping (no single condition of interest) | In DE, PCs can absorb the condition; SVA and RUV protect it |

  * *Why qSVA in postmortem brain:* when RNA quality correlates with diagnosis, quality metrics alone do not remove the
    confounding; qSVA, based on degradation-associated regions, more than tripled replication across independent
    studies [Jaffe 2017].
  * *Why PCA only for QTL:* for molecular QTL mapping PCA matched or outperformed SVA, PEER (probabilistic estimation of expression residuals) and HCP (hidden covariates with prior) and is much faster
    [Zhou 2022; Stegle 2012], but in DE the factors must be estimated with the condition protected.
* **Number of factors**: use the method's estimate (`num.sv` in SVA; for RUV, choose $k$ with relative log expression
  plots and p-value histograms), and report how the DE list changes over a range of $k$.
  * *Why:* too many factors absorb real signal; too few leave confounding. Stability over $k$ shows whether
    conclusions depend on the choice.
  * Recipe: `recipes/bulk/hidden_factors.py` (principal components of residuals with the known design protected;
    DE across a grid of $k$ with overlap against the chosen $k$).
* MUST NOT use batch-corrected matrices for testing [Nygaard 2015†; Zhang 2020†].

## 8. Differential expression

* **Default in Python: `pydeseq2`** [Muzellec 2023]; R options: DESeq2 [Love 2014†; Love 2016a], edgeR quasi-likelihood
  [Chen 2016], limma-voom [Law 2014†]. When `pydeseq2` results are reported, verify on one contrast that they agree with
  DESeq2.
* **Wald test** for a two-level contrast; **likelihood-ratio test** for factors with several levels or time courses.
  `pydeseq2` 0.5 has no likelihood-ratio test; run it with DESeq2 in R through a script (`CODING_GUIDELINES.md` §3,
  languages and tooling).
* **Shrunken log fold changes** for ranking and plots [Zhu 2018†]; keep unshrunken estimates and standard errors in the
  results table.
* **Minimum effect**: test against a threshold (`lfcThreshold`, TREAT) rather than filtering by fold change afterwards
  [McCarthy 2009†; Ebrahimpoor 2021†].
* **Large human cohorts**: check FDR calibration by repeating the analysis with permuted condition labels
  [Li 2022†].
  * Recipe: `recipes/bulk/fdr_calibration.py`.
* Results table: `gene_id`, `gene_symbol`, `base_mean`, `log2_fold_change`, `lfc_se`, `log2_fold_change_shrunk`,
  `stat`, `p_value`, `p_adjusted`, contrast name.
* Recipe: `recipes/bulk/deseq2.py` (pydeseq2 with this results table, minimum-effect test, shrinkage and VST).

## 9. Transcript-level analyses

* Distinguish DGE (differential gene expression), DTE (differential transcript expression) and DTU (differential
  transcript usage: a change in isoform proportions).
* **DTE and DTU use inferential replicates** (§4): swish [Zhu 2019] or edgeR/limma with divided counts [Baldoni 2023;
  Baldoni 2025]. Without replicates, DRIMSeq or DEXSeq with stage-wise testing [Love 2018; Anders 2012] and isoform
  prefiltering [Soneson 2016b].
  * *Why:* methods that ignore quantification uncertainty inflate the FDR for transcripts with high ambiguity
    [Zhu 2019].

## 10. Enrichment and activity inference

* **Ranked tests (GSEA, gene set enrichment analysis)** on all tested genes ranked by a signed statistic [Subramanian 2005; Korotkevich 2016;
  Fang 2022].
* **Over-representation analysis** with the tested genes as background (`STATISTICS_GUIDELINES.md` §9)
  [Wijesooriya 2022†].
* **Length and detection bias**:
  * goseq adjusts over-representation analysis for the higher detection power of long genes [Young 2010]; it applies
    to full-length protocols, not to 3' tag data;
  * sample-specific length bias also affects ranked tests (§5) [Mandelboum 2019].
* **Inter-gene correlation**: tests that permute genes assume genes are independent; correlated genes in a set inflate
  type I error. camera estimates the correlation and keeps the error controlled [Wu 2012]; under phenotype permutation,
  several methods showed large false-positive rates [Tarca 2013]. Methods differ in runtime and in which terms they
  prioritize [Geistlinger 2020].
* **Pathway and transcription-factor activity ("footprints")** with `decoupler` [Badia-i-Mompel 2022]:
  * pathways with PROGENy, which infers activity from downstream perturbation-response genes rather than pathway
    members [Schubert 2018];
  * transcription factors with CollecTRI regulons, which outperformed other resources [Müller-Dott 2023]; curated
    interactions perform best [Garcia-Alonso 2019];
  * the univariate linear model (ULM) as default method; the gene-set resource matters more than the method
    [Holland 2020].
* Record the gene-set database and release; for mouse, use mouse collections or one-to-one orthologs.
* A practical protocol for the full workflow is given in [Reimand 2019].

## 11. Human and mouse

* Species-specific annotation; never mix human and mouse IDs.
* Cross-species comparisons use one-to-one orthologs from a stated release; genes without them are reported
  separately.

## 12. What to report

In addition to `STATISTICS_GUIDELINES.md` §12:

1. Genome build, annotation release, pipeline and tool versions, Salmon options.
2. Per-sample QC table, sex and identity check results.
3. tximport mode, filtering rule, normalization, design formula, contrasts, known and hidden factors with the method
   and number.
4. DE tool, test, threshold test, shrinkage.
5. Enrichment and activity methods, ranking statistic or background, resources and releases.

## 13. Implementation

| Task | Python | R reference (via script) |
|---|---|---|
| FASTQ to counts | `nf-core/rnaseq` (Nextflow) | — |
| QC aggregation | MultiQC | — |
| Transcript-to-gene | — | `tximport` |
| Differential expression, VST | `pydeseq2` | DESeq2, edgeR, limma |
| Known/hidden factors | — | `sva`, `RUVSeq`, `ruv`, qSVA (`qsvaR`) |
| Length bias | — | `cqn`, `EDASeq` |
| Transcript expression and usage | — | `fishpond` (swish), edgeR (`catchSalmon`), DRIMSeq, DEXSeq, stageR |
| GSEA / ORA (over-representation analysis) | `gseapy` | `fgsea`, `goseq`, `limma::camera` |
| Pathway and transcription factor (TF) activity | `decoupler` (PROGENy, CollecTRI) | `decoupleR` |
