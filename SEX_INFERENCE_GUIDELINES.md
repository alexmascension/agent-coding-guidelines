# Sex inference guidelines

> Source: <https://github.com/alexmascension/agent-coding-guidelines> · Version 2.0.0 (2026-10-15)
>
> Key words **MUST**, **SHOULD**, **MAY** follow RFC 2119. Bracketed keys such as `[Olney 2020]` point to
> `OMICS_REFERENCES.md`; keys with a dagger point to `STATISTICS_REFERENCES.md`. Each rule is followed by its
> rationale: ***Why*** (mechanism), ***Trade-off*** (cost) and ***Deviate when*** (justified exceptions), where they
> apply.
> A rule without a key word is a default with the force of SHOULD. Rules carry stable identifiers (e.g.
> `SEX-WHY-001`); overrides refer to them (`CODING_GUIDELINES.md` §0, rule identifiers).

## 0. Scope and purpose

* Infer the sex-chromosome complement of each sample from transcriptomic data, compare it with the reported sex, and
  use it to (a) detect sample swaps and annotation errors, (b) choose sex-matched references, (c) flag anomalies
  (aneuploidy, mosaic loss of Y, mixtures).
* Applies to bulk RNA-seq and to single-cell and single-nucleus RNA-seq, human and mouse. The procedure is independent
  of aligner and counting tool; only the inputs (counts per marker gene and library size) are required.
* Copy with `OMICS_REFERENCES.md`. Related: `SAMPLE_IDENTITY_GUIDELINES.md` (genotype-based identity),
  `RNASEQ_BULK_GUIDELINES.md` §2 (reference genome).
* Terminology: this guide infers **chromosomal sex** (XX, XY and anomalies) from expression; it does not infer gender.
  Results are labelled `XX`, `XY`, etc., never "female" or "male", and are compared with the reported sex in a separate
  column.

## 1. Why run it on every dataset

* `SEX-WHY-001` **Every dataset with sex metadata MUST be checked**, including public datasets and the user's own.
  * *Why:* when reported sex was compared with sex-specific gene expression, 46% of public human transcriptomics
    datasets contained at least one apparently mislabelled sample, and in a single cohort two of four datasets had
    mix-ups, pointing to laboratory errors rather than typos [Toker 2016]. A sex mismatch is the cheapest detectable
    symptom of a sample swap; genotype checks find more (`SAMPLE_IDENTITY_GUIDELINES.md`).
  * *Trade-off:* minutes of computation per dataset.
  * *Scope:* no exceptions for analyses that use sex as a variable or that combine samples across donors.

## 2. Marker genes

* `SEX-MARK-001` **Female-specific signal: *XIST*** (mouse *Xist*), optionally with *TSIX*.
  * *Why:* *XIST* is expressed from the inactive X, so in normal adult tissues it is present in XX cells and
    essentially absent in XY cells. It is nuclear-retained, which makes it especially strong in single-nucleus data.
  * *Limits:* *XIST* is also expressed with more than one X (e.g. XXY) and in some tumours of male patients (e.g.
    testicular germ-cell tumours); *XIST* in a sample with Y signal is a reason for review (`XIST+Y+`, §5), not
    proof of a swap.
  * *Deviate when:* do not add other X-linked genes as female markers. No other adult X gene is XX-specific; genes
    that escape X inactivation (*KDM5C*, *KDM6A*, *DDX3X*, *ZFX*, *EIF1AX*) are only modestly higher in XX and overlap
    between sexes [Tukiainen 2017; Oliva 2020].
* `SEX-MARK-002` **Male-specific signal: broadly expressed male-specific region of Y (MSY) genes**, summed. Human: *RPS4Y1*,
  *DDX3Y*, *KDM5D*, *UTY*, *EIF1AY*, *USP9Y*, *ZFY* (and *NLGN4Y* in brain). Mouse: *Ddx3y*, *Eif2s3y*, *Kdm5d*, *Uty*.
  * *Why:* a five-gene panel (*RPS4Y1*, *EIF1AY*, *DDX3Y*, *KDM5D*, *XIST*) classified sex across 24 tissue collections
    and platforms [Staedtler 2013]; X-degenerate Y genes are expressed across most tissues [Godfrey 2020]. Summing
    several genes reduces dependence on any one gene's tissue-specific expression.
  * *Deviate when:* drop genes not expressed in the tissue at hand (check on samples of known sex).
* `SEX-MARK-003` **Within-sample ratio from X/Y homolog pairs**: for each pair, $r = Y/(X+Y)$, then the median across pairs. Human
  pairs include *RPS4X/RPS4Y1*, *DDX3X/DDX3Y*, *KDM5C/KDM5D*, *KDM6A/UTY*, *EIF1AX/EIF1AY*, *USP9X/USP9Y*, *ZFX/ZFY*;
  mouse: *Ddx3x/Ddx3y*, *Eif2s3x/Eif2s3y*, *Kdm5c/Kdm5d*, *Kdm6a/Uty*.
  * *Why:* the X homolog is a built-in denominator. The ratio cancels factors that scale both homologs equally
    (library size, overall depth) and separates XY ($r$ substantial) from XX ($r \approx 0$) without a cohort-level
    threshold. It is also the most sensitive signal for male material mixed into an XX sample.
  * *Limits:* tissue, RNA quality and mapping can affect the two homologs differently, so the ratio is more robust
    than absolute values but not invariant to them; inspect its distribution per tissue and batch.
  * *Trade-off:* homologs share sequence; with a reference that does not mask the Y PAR (pseudoautosomal region), or with multi-mapping reads,
    reads can be misassigned (§4).

## 3. Counting and normalization

* `SEX-COUNT-001` **Count marker genes over the whole gene body (exons and introns)** with uniquely mapped reads, when the library
  captures intronic reads (total RNA, nuclei).
  * *Why:* in total RNA and in nuclei about half of the reads are intronic, and several Y genes have long loci, so the
    signal multiplies several-fold. Introns of X-degenerate Y genes have diverged from their X homologs, so unique
    reads are unambiguous.
  * *Deviate when:* for poly(A) or 3'-tag libraries, exonic counts are sufficient.
* `SEX-COUNT-002` **Normalize to counts per million uniquely mapped reads (or fragments)** and analyse on the $\log_{10}(1+x)$ scale.
  * *Why:* the distributions are bimodal with many values near zero; the log scale separates the modes and stabilizes
    Poisson noise.
* `SEX-COUNT-003` **Screening depth:** for XX versus XY, about 1 million read pairs per sample is a starting value. Expected values are on the
  order of tens to hundreds of informative reads, and observing zero MSY reads in an XY sample is then negligible.
  Re-evaluate anomalies at full depth.
  * *Why:* detecting a 5–10% mixture or partial loss of Y needs the full library; a 10% male contribution to an XX
    sample adds only a few Y reads at 1 million pairs.
* `SEX-COUNT-004` **Minimum informative reads:** samples with fewer than a fixed number of informative reads (XIST + MSY + X homologs;
  about 30 as a starting value) are `LOW_SIGNAL` and are re-run at higher depth, not called.

## 4. Reference genome

* `SEX-REF-001` **Routing rule (shared with `RNASEQ_BULK_GUIDELINES.md` §2):** sex-chromosome-complement-informed references SHOULD
  be used for the main alignment whenever sex-chromosome or PAR genes are analysed or sex is a variable of the study,
  and MAY be used otherwise. Y is fully masked for samples called `XX`; the Y PAR is masked for all other samples (the
  PAR is then represented only on X) [Olney 2020]. Uncertain calls (`XY_lowY`, `XIST+Y+`, `XIST-Y-`, `LOW_SIGNAL`)
  go to the XY reference, which keeps any Y reads, and stay flagged.
  * *Why:* with the default reference, reads from X/Y homologous and PAR sequence are split or misassigned, which
    lowers PAR expression estimates and produces spurious Y expression in XX samples [Olney 2020; Webster 2019].
  * *Trade-off:* two references and a routing step based on the screening call (§6).
* `SEX-REF-002` The screening itself MUST use a reference that contains Y (XY reference); otherwise Y signal cannot be observed.
* `SEX-REF-003` The annotation MUST NOT contain features on masked sequence (e.g. genes annotated in the Y PAR for the XY reference);
  remove them, or the corresponding transcripts are all `N`.

## 5. Thresholds and categories

* `SEX-THR-001` **Thresholds are data-driven by default**: on the log scale, a two-group split (Otsu, or a two-component Gaussian
  mixture), with the cut placed between the medians of the two groups. Unsupervised clustering of Y-gene expression is
  the same idea applied to microarrays [Buckberry 2014].
  * *Why:* absolute levels depend on tissue, protocol and depth, so a fixed threshold that works in blood fails in
    brain.
  * *Trade-off / Deviate when:* a data-driven split needs both groups present in reasonable numbers. With a single-sex
    cohort, a strongly unbalanced cohort, or fewer than about 20 samples, use thresholds calibrated on samples of known
    sex processed the same way, and say so.
* `SEX-THR-002` **Y presence requires both** absolute Y signal above threshold **and** a median homolog ratio above a minimum (e.g.
  $r \ge 0.10$).
  * *Why:* requiring both avoids calling Y from mapping noise in XX samples.
* `SEX-THR-003` **Categories** (one per sample):

  | Call | XIST | Y | Meaning | Action |
  |---|---|---|---|---|
  | `XX` | high | absent | XX | XX reference |
  | `XY` | absent | present | XY | XY reference |
  | `XY_lowY` | absent | present but below a fraction (e.g. 0.5) of the XY median | Possible mosaic loss of Y, or mixture | XY reference; review (§7) |
  | `XIST+Y+` | high | present | XXY, or mixture of XX and XY material | XY reference; review with genotype |
  | `XIST-Y-` | absent | absent | X0, or technical failure | XY reference; review |
  | `LOW_SIGNAL` | — | — | Too few informative reads | Re-run deeper |

* `SEX-THR-004` **Flag `possible_mix`** in `XX` samples when the homolog ratio is above a trace level (e.g. 0.02) or Y signal is far
  above the XX median.
* `SEX-THR-005` **Report, for every call, the distance to each threshold** (in log units), not only the category.
  * *Why:* calls near a threshold are the ones that need review; a bare label hides them.
* `SEX-THR-006` **Calibrate the mixture and loss-of-Y flags with in-silico mixtures**: combine reads from an XX and an XY sample at
  known fractions (e.g. 1, 2, 5, 10, 20%) and record the fraction at which each flag fires.
  * *Why:* the depth, minimum-read and ratio values of this guide (1 million pairs, $r \ge 0.10$, 0.02, 0.5 of the
    XY median) are starting values, not guaranteed sensitivities; their real behaviour depends on depth, tissue and
    library type. They MUST be calibrated (in-silico mixtures, samples of known karyotype) before a project relies on
    them, and the calibrated values are recorded.

## 6. Combining evidence and deciding

* `SEX-DEC-001` **Order of evidence** for the karyotype used in analysis and for reference routing:
  1. WGS (whole-genome sequencing) or genotyping-array karyotype of the donor, when available.
  2. Expression call (this guide).
  3. Reported sex.
* `SEX-DEC-002` **Genotype-based sex** from the same RNA BAM: heterozygosity on non-PAR X (XX and XXY have heterozygous sites; XY
  does not) and Y depth, e.g. with somalier [Pedersen 2020] or peddy (DNA) [Pedersen 2017].
  * *Why:* it is orthogonal to expression; a sample that is `XX` by expression but has no X heterozygosity is
    suspicious.
  * *Trade-off:* in RNA, only expressed genes are covered, and allele-specific expression (including X inactivation in
    clonal tissue) can hide heterozygosity; treat it as supporting evidence.
* `SEX-DEC-003` **Expression calls are screening results, not karyotypes.** Keep separate columns: reported sex (metadata),
  expression call (this guide), chromosome complement confirmed by DNA (only from WGS, arrays or DNA-based
  genotype), flags (`possible_mix`, `XY_lowY`, review status) and the reference used. A complement is "confirmed" only
  with DNA evidence.
* `SEX-DEC-004` **All discordances are flagged; they MUST NOT be corrected silently.** Metadata changes are recorded in the decisions log with
  the evidence.
* `SEX-DEC-005` **Discordance between samples of one donor MUST trigger review; it MUST NOT by itself trigger reassignment or a
  metadata change.** A sample swap is the most frequent explanation (`SAMPLE_IDENTITY_GUIDELINES.md`), but mosaic or
  tissue-specific loss of Y (§7), tumour or pathology, low quality, and sex-mismatched transplantation (e.g. bone
  marrow, whose cells reach blood and tissue-resident myeloid populations) also produce discordant calls. Decide with
  genotype identity and DNA evidence, and record the decision.

## 7. Loss of Y and other anomalies

* `SEX-LOY-001` **Mosaic loss of chromosome Y (LOY) is common and biological.** About 20% of men in UK Biobank have detectable LOY in
  blood [Thompson 2019], and it is associated with mortality [Forsberg 2014] and disease [Sano 2022]. In brain, LOY is
  enriched in microglia and more frequent in Alzheimer's disease [Vermeulen 2022]. It can be scored from bulk and
  single-cell RNA-seq [Dumanski 2021].
  * *Consequence:* `XY_lowY` in aged tissue MUST NOT be interpreted as contamination by default. Interpret it with cell
    composition (e.g. microglial fraction from deconvolution) and, when single-cell data exist, with per-cell Y
    detection.
* `SEX-LOY-002` `XIST+Y+` is resolved with genotype: XXY shows X heterozygosity and a single genotype; a mixture shows excess
  heterozygosity, allele fractions inconsistent with one individual, and high relatedness to another sample.
* `SEX-LOY-003` `XIST-Y-` with good library quality and no X heterozygosity suggests X0; with poor quality, a failed library.

## 8. Single-cell and single-nucleus data

* `SEX-SC-001` **Call sex per sample (per donor), from pseudobulk counts**, using §2–5.
  * *Why:* per-cell detection is sparse (many cells have zero counts of every marker), so per-cell calls are noisy;
    the sample-level signal is robust.
* `SEX-SC-002` **Per-cell scores MAY be used** to (a) check that a sample is homogeneous, (b) detect cells from another donor in
  pooled designs, (c) estimate per-cell LOY in XY donors.
* `SEX-SC-003` **Ambient RNA carries marker transcripts**: in pools of XX and XY donors, ambient *XIST* and Y transcripts appear in
  every droplet. Use ambient-corrected counts and, for pools, genotype demultiplexing rather than sex markers
  (`SAMPLE_IDENTITY_GUIDELINES.md` §6).
* `SEX-SC-004` In nuclei, *XIST* is strong (nuclear RNA) and Y genes are often detected from intronic reads; count with introns.

## 9. Mouse

* `SEX-MOUSE-001` Use *Xist* and the mouse MSY genes and homolog pairs of §2. Mouse Y carries few broadly expressed genes, so the Y sum
  is smaller than in human.
* `SEX-MOUSE-002` Inbred strains have little heterozygosity, so genotype-based sex from X heterozygosity does not work; rely on
  expression (and on breeding records).
* `SEX-MOUSE-003` Sex-chromosome complement affects autosomal gene expression independently of gonadal sex [Wijchers 2010]; models
  such as XY females are not identified by these markers alone.

## 10. Outputs and report

* `SEX-REP-001` One table per dataset with, per sample: reported sex, call, distance to each threshold, informative reads,
  normalized XIST and Y, median homolog ratio, `possible_mix`, chosen reference, concordance (`concordant`,
  `DISCORDANT`, `no_reported_sex`, `uncertain_call`), and review flag.
* `SEX-REP-002` A scatter plot of $\log_{10}(1+\text{XIST})$ against $\log_{10}(1+\text{Y})$ coloured by reported sex, with the
  thresholds, discordant samples circled, and a histogram of the homolog ratio.
* `SEX-REP-003` The thresholds used (and whether data-driven or calibrated), the marker lists and the reference versions.
* `SEX-REP-004` Counts of discordances and anomalies, and the actions taken for each.
